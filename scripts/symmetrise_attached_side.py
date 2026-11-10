#!/usr/bin/env python3
"""Mirror-symmetric start that COPIES one half onto the other (instead of averaging the two halves).

The half kept is the more attached one: higher peak electron temperature at its two targets (the target columns
are the first and last physical column of each of the two reflection blocks), reported together with the peak
target density; override with side='upper'|'lower'.  Plasma state (b2fstate) and the short virtual-neutral
background (fort.13) are treated; parallel velocity (odd) changes sign, poloidal face fluxes follow the same
face-reflection rule as symmetrise_state.py.  Species inventories and energies are NOT conserved (the other
half's content is replaced); the report gives the before/after totals.

usage: symmetrise_attached_side.py <b2fstate> <b2fgmtry> <fort.13> <baserun dir> <out b2fstati> <out fort.13> [auto|upper|lower]
"""
import json, math, sys
from array import array
from pathlib import Path
import numpy as np
from scipy.io import FortranFile
from audit_mesh_symmetry import cf_fields
from symmetrise_state import read_state, write_state


def project_state(source, geometry, destination, side='auto', proton_mass=1.67262192369e-27):
    version, fields = read_state(source)
    data = {name: values for kind, count, name, values in fields}
    nx, ny, ns = data['nx,ny,ns']
    n = (nx+2)*(ny+2)
    geo = cf_fields(geometry, {'nx,ny', 'leftix', 'rightix', 'rightiy', 'cry', 'vol'})
    cuts = [x for x in range(-1, nx+1) if geo['leftix'][x+1+(nx+2)] == -2]
    if len(cuts) != 2 or cuts[0] != -1:
        raise ValueError('Expected two reflection blocks; found '+str(cuts))
    blocks = [(-1, cuts[1]-1), (cuts[1], nx)]
    xmap = {x: lo+hi-x for lo, hi in blocks for x in range(lo, hi+1)}
    mapping = []
    for y in range(-1, ny+1):
        for x in range(-1, nx+1):
            mapping.append(xmap[x]+1+(nx+2)*(y+1))
    zc = [sum(geo['cry'][i+n*k] for k in range(4))/4 for i in range(n)]
    upper = [z > 0 for z in zc]
    pairs = [(i, j) for i, j in enumerate(mapping) if i < j]
    # attachment metric per half: sum over its two target columns of the peak Te (eV)
    ev = 1.602176634e-19
    tcols = [0, cuts[1]-2, cuts[1]+1, nx-1]
    score = {True: [], False: []}
    dens = {True: [], False: []}
    for x in tcols:
        for half in (True, False):
            cells = [x+1+(nx+2)*(y+1) for y in range(ny) if upper[x+1+(nx+2)*(y+1)] == half]
            if cells:
                score[half].append(max(data['te'][i] for i in cells)/ev)
                dens[half].append(max(data['ne'][i] for i in cells))
    metric = {'upper_peak_Te_eV_by_target_column': score[True], 'lower_peak_Te_eV_by_target_column': score[False],
              'upper_peak_ne_by_target_column': dens[True], 'lower_peak_ne_by_target_column': dens[False]}
    if side == 'auto':
        keep_upper = sum(score[True]) >= sum(score[False])
    else:
        keep_upper = side == 'upper'
    original = {k: array('d', data[k]) for k in ['na', 'ne', 'ua', 'te', 'ti', 'kinrgy', 'po']}
    ov = geo['vol']
    physical = [x+1+(nx+2)*(y+1) for y in range(ny) for x in range(nx) if x not in [cuts[1]-1, cuts[1]]]
    charge = [.5*(a+b) for a, b in zip(data['zamin'], data['zamax'])]
    ions = [s for s, z in enumerate(charge) if z > 0]

    def totals(state):
        inv = [math.fsum(state['na'][i+n*s]*ov[i] for i in physical) for s in range(ns)]
        ee = math.fsum(1.5*state['ne'][i]*state['te'][i]*ov[i] for i in physical)
        return inv, ee
    before = totals(original)
    kept, other = {}, {}
    for i, j in pairs:
        c, o = (i, j) if upper[i] == keep_upper else (j, i)
        if upper[c] != keep_upper or upper[o] == keep_upper:
            raise ValueError('Mirror pair not split between the halves: '+str((i, j)))
        kept[o] = c
        for s in range(ns):
            data['na'][o+n*s] = original['na'][c+n*s]
            data['ua'][o+n*s] = -original['ua'][c+n*s]
            data['kinrgy'][o+n*s] = original['kinrgy'][c+n*s]
        for k in ['ne', 'te', 'ti', 'po']:
            data[k][o] = original[k][c]
    handled = {'na', 'ne', 'ua', 'te', 'ti', 'po', 'kinrgy'}
    metadata = {'nx,ny,ns', 'label', 'zamin', 'zamax', 'zn', 'am', 'time'}
    face_fields = []
    for kind, count, name, values in fields:
        if name in handled or name in metadata:
            continue
        old = array('d', values)
        if name == 'uadia':
            for block in range(2*ns):
                sign = -1 if block % 2 == 0 else 1
                for o, c in kept.items():
                    values[o+n*block] = sign*old[c+n*block]
            continue
        if name == 'fch_p':
            blocks_n = 1
        elif count % (2*n) == 0 and (name.startswith('f') or name in ['vaecrb', 'vadia', 'wadia', 'veecrb', 'vedia']):
            blocks_n = count//n
        else:
            raise ValueError('Unclassified restart field: '+name)
        face_fields.append(name)
        for block in range(blocks_n):
            poloidal = name == 'fch_p' or block % 2 == 0
            sign = -1 if poloidal else 1
            partner = []
            for i, cell_partner in enumerate(mapping):
                if poloidal:
                    if geo['leftix'][i] == -2:
                        partner.append(None)
                        continue
                    x, y = geo['rightix'][cell_partner], geo['rightiy'][cell_partner]
                    if not (-1 <= x <= nx and -1 <= y <= ny):
                        raise ValueError('Reflected physical face is outside geometry')
                    partner.append(x+1+(nx+2)*(y+1))
                else:
                    partner.append(cell_partner)
            done = set()
            for i, j in enumerate(partner):
                if j is None:
                    values[i+n*block] = 0.
                    continue
                if i in done:
                    continue
                if j == i:
                    values[i+n*block] = 0. if poloidal else old[i+n*block]
                    done.add(i)
                    continue
                ki, kj = upper[i] == keep_upper, upper[j] == keep_upper
                if ki and not kj:
                    values[j+n*block] = sign*old[i+n*block]
                elif kj and not ki:
                    values[i+n*block] = sign*old[j+n*block]
                else:  # ambiguous owner (face across the mirror plane): average, as symmetrise_state.py
                    v = .5*(old[i+n*block]+sign*old[j+n*block])
                    values[i+n*block], values[j+n*block] = v, sign*v
                done.add(i)
                done.add(j)
    for name in ['na', 'ne', 'te', 'ti']:
        if min(data[name]) <= 0:
            raise ValueError('Nonpositive '+name)
    after = totals(data)
    write_state(destination, version, fields)
    report = {'kept_half': 'upper' if keep_upper else 'lower', 'rule': 'auto: larger summed peak target Te' if side == 'auto' else side,
              'attachment_metric': metric, 'inventory_relative_change': [abs(a-b)/abs(a) for a, b in zip(before[0], after[0])],
              'electron_thermal_energy_relative_change': abs(before[1]-after[1])/before[1],
              'cache_face_fields_projected': face_fields, 'input': str(source), 'output': str(destination)}
    destination.with_name('state_projection_attached.json').write_text(json.dumps(report, indent=2)+'\n')
    return report, keep_upper


def project_neutrals(source, base, destination, keep_upper):
    """Copy the kept half of the short virtual-neutral record (9 species, 5 moments per triangle) with parity."""
    mapping = np.loadtxt(base/'triangle_reflection.txt', dtype=int)
    tok = (base/'fort.33').read_text().split()
    nv = int(tok[0]); r = np.array(tok[1:1+nv], float); z = np.array(tok[1+nv:1+2*nv], float)
    tri = (base/'fort.34').read_text().split()
    nt = int(tri[0]); v = np.array(tri[1:1+4*nt], int).reshape(nt, 4)[:, 1:]-1
    zt = z[v].mean(1)
    with FortranFile(source, 'r') as f:
        a = f.read_reals(np.float64)
    nc = a.size//45
    assert nc == len(mapping)+6 and nt == len(mapping), (nc, len(mapping), nt)
    m = np.r_[mapping, np.arange(len(mapping), nc)]
    assert np.array_equal(m[m], np.arange(nc))
    old = a.reshape((9, nc, 5), order='F'); new = old.copy()
    kept = 0
    for i in range(len(mapping)):
        j = m[i]
        if i < j:
            c, o = (i, j) if (zt[i] > 0) == keep_upper else (j, i)
            if (zt[c] > 0) != keep_upper or (zt[o] > 0) == keep_upper:
                continue  # pair straddles Z=0 or is degenerate: leave as it is
            new[:, o, :] = old[:, c, :]
            new[:, o, 3:] *= -1
            kept += 1
    with FortranFile(destination, 'w') as f:
        f.write_record(new.ravel(order='F'))
    rep = {'pairs_copied': kept, 'pairs_total': int((np.arange(len(mapping)) < mapping).sum()), 'source': str(source),
           'destination': str(destination), 'cells': nc, 'kept_half': 'upper' if keep_upper else 'lower'}
    destination.with_name('neutral_projection_attached.json').write_text(json.dumps(rep, indent=2)+'\n')
    return rep


if __name__ == '__main__':
    st, geo, fort13, base, out_state, out_n = map(Path, sys.argv[1:7])
    side = sys.argv[7] if len(sys.argv) > 7 else 'auto'
    rep, ku = project_state(st, geo, out_state, side)
    print(json.dumps(rep, indent=2))
    print(json.dumps(project_neutrals(fort13, base, out_n, ku), indent=2))
