#!/usr/bin/env python3
"""Project a STEP restart onto reflection symmetry using conserved cell totals.

Preserves each species inventory, electron thermal energy, and combined charged
ion thermal/parallel kinetic energy over physical cells. Removes the forbidden
even component of parallel momentum; reflects cached face fluxes with their
proper face index shift. This is preparation, not a converged solution.
"""
from array import array
import argparse
import json
import math
from pathlib import Path
import re

from audit_mesh_symmetry import cf_fields


def read_state(path):
    header = re.compile(r'^\*cf:\s+(int|real|char)\s+(\d+)\s+(\S+)')
    fields = []
    with path.open() as stream:
        version = next(stream)
        current, lines = None, []

        def finish():
            if current is None:
                return
            kind, count, name = current
            if kind == 'char':
                values = ''.join(lines)
            else:
                values = array('q' if kind == 'int' else 'd')
                for line in lines:
                    values.extend((int(v) if kind == 'int' else float(v.replace('D', 'E')))
                                  for v in line.split())
                if len(values) != count:
                    raise ValueError('Wrong field length: '+name)
            fields.append([kind, count, name, values])

        for line in stream:
            match = header.match(line)
            if match:
                finish()
                kind, count, name = match.groups()
                current, lines = (kind, int(count), name), []
            else:
                lines.append(line)
        finish()
    if len({row[2] for row in fields}) != len(fields):
        raise ValueError('Duplicate state fields')
    return version, fields


def write_state(path, version, fields):
    """Match the native cfruin/cfrure layout (6 values of 22 characters, 13 decimals), so that a b2fstate written by b2mn
    is never smaller than this file (the run scripts judge a run by file size)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as stream:
        stream.write(version)
        for kind, count, name, values in fields:
            stream.write(f'{"*cf:":8s}{kind:8s}{count:12d}    {name:32s}\n')
            if kind == 'char':
                stream.write(values)
            elif kind == 'int':
                for i in range(0,count,12):
                    stream.write(''.join(f' {v:10d}' for v in values[i:i+12])+'\n')
            else:
                for i in range(0,count,6):
                    formatted=[]
                    for value in values[i:i+6]:
                        mantissa, exponent=f'{value:.13E}'.split('E')
                        literal=f'{mantissa}E{int(exponent):+04d}'
                        formatted.append(f'{literal:>22s}')
                    stream.write(''.join(formatted)+'\n')


def project(source, old_geometry, new_geometry, destination, proton_mass=1.67262192369e-27):
    version, fields = read_state(source)
    data = {name: values for kind, count, name, values in fields}
    nx, ny, ns = data['nx,ny,ns']
    n = (nx+2)*(ny+2)
    keys = {'nx,ny', 'vol', 'rightix', 'rightiy', 'leftix'}
    old_geo, geo = cf_fields(old_geometry, keys), cf_fields(new_geometry, keys)
    if old_geo['nx,ny'] != [nx, ny] or geo['nx,ny'] != [nx, ny]:
        raise ValueError('Geometry/state shapes disagree')
    # Two blocks, each reflected onto itself; the second starts at the second column without a west neighbour (iy=0).
    cuts = [x for x in range(-1, nx+1) if geo['leftix'][x+1+(nx+2)] == -2]
    if len(cuts) != 2 or cuts[0] != -1:
        raise ValueError('Expected two reflection blocks; found block starts '+str(cuts))
    blocks = [(-1, cuts[1]-1), (cuts[1], nx)]
    xmap = {x: lo+hi-x for lo, hi in blocks for x in range(lo, hi+1)}
    pairs, mapping = [], []
    for y in range(-1, ny+1):
        for x in range(-1, nx+1):
            i, j = x+1+(nx+2)*(y+1), xmap[x]+1+(nx+2)*(y+1)
            mapping.append(j)
            if i < j:
                pairs.append((i, j))
    physical = [x+1+(nx+2)*(y+1) for y in range(ny) for x in range(nx) if x not in [cuts[1]-1, cuts[1]]]
    ov, nv = old_geo['vol'], geo['vol']
    charge = [.5*(a+b) for a, b in zip(data['zamin'], data['zamax'])]
    ions = [s for s, z in enumerate(charge) if z > 0]
    # EFGW Fundamental-Constants/0.1.1 supplies CODATA 2018 for this build.
    # Pass the matching build constant when using another environment.
    original = {k: array('d', data[k]) for k in ['na', 'ne', 'ua', 'te', 'ti', 'kinrgy']}

    def totals(state, volumes):
        inventories = [math.fsum(state['na'][i+n*s]*volumes[i] for i in physical) for s in range(ns)]
        electrons = math.fsum(1.5*state['ne'][i]*state['te'][i]*volumes[i] for i in physical)
        ion_energy = math.fsum(volumes[i]*math.fsum(state['na'][i+n*s]*
            (1.5*state['ti'][i]+state['kinrgy'][i+n*s]) for s in ions) for i in physical)
        return inventories, electrons, ion_energy

    before = totals(original, ov)
    for i, j in pairs:
        denominator = nv[i]+nv[j]
        if min(ov[i], ov[j], nv[i], nv[j]) <= 0:
            raise ValueError('Nonpositive cell volume')
        for s in range(ns):
            a, b = i+n*s, j+n*s
            particles = original['na'][a]*ov[i]+original['na'][b]*ov[j]
            density = particles/denominator
            velocity = (original['na'][a]*ov[i]*original['ua'][a]-
                        original['na'][b]*ov[j]*original['ua'][b])/particles
            data['na'][a] = data['na'][b] = density
            data['ua'][a], data['ua'][b] = velocity, -velocity
            data['kinrgy'][a] = data['kinrgy'][b] = .5*data['am'][s]*proton_mass*velocity**2
        density_e = math.fsum(charge[s]*data['na'][i+n*s] for s in ions)
        data['ne'][i] = data['ne'][j] = density_e
        electron_energy = (original['ne'][i]*original['te'][i]*ov[i]+
                           original['ne'][j]*original['te'][j]*ov[j])
        data['te'][i] = data['te'][j] = electron_energy/(density_e*denominator)
        old_ion_energy = math.fsum(original['na'][cell+n*s]*ov[cell]*
            (1.5*original['ti'][cell]+original['kinrgy'][cell+n*s]) for cell in [i, j] for s in ions)
        new_kinetic = denominator*math.fsum(data['na'][i+n*s]*data['kinrgy'][i+n*s] for s in ions)
        ion_density = math.fsum(data['na'][i+n*s] for s in ions)
        temperature = (old_ion_energy-new_kinetic)/(1.5*ion_density*denominator)
        if temperature <= 0:
            raise ValueError('Energy projection produced nonpositive ion temperature')
        data['ti'][i] = data['ti'][j] = temperature
        data['po'][i] = data['po'][j] = .5*(data['po'][i]+data['po'][j])

    handled = {'na', 'ne', 'ua', 'te', 'ti', 'po', 'kinrgy'}
    metadata = {'nx,ny,ns', 'label', 'zamin', 'zamax', 'zn', 'am', 'time'}
    face_fields = []
    for kind, count, name, values in fields:
        if name in handled or name in metadata:
            continue
        old = array('d', values)
        if name == 'uadia':
            if count != 2*n*ns:
                raise ValueError('Unexpected uadia shape')
            for block in range(2*ns):
                sign = -1 if block % 2 == 0 else 1
                for i, j in pairs:
                    value = .5*(old[i+n*block]+sign*old[j+n*block])
                    values[i+n*block], values[j+n*block] = value, sign*value
            continue
        if name == 'fch_p':
            if count != n:
                raise ValueError('Unexpected fch_p shape')
            blocks = 1
        elif count % (2*n) == 0 and (name.startswith('f') or name in ['vaecrb', 'vadia', 'wadia', 'veecrb', 'vedia']):
            blocks = count//n
        else:
            raise ValueError('Unclassified restart field: '+name)
        face_fields.append(name)
        for block in range(blocks):
            poloidal = name == 'fch_p' or block % 2 == 0
            sign = -1 if poloidal else 1
            for i, cell_partner in enumerate(mapping):
                if poloidal:
                    if geo['leftix'][i] == -2:
                        values[i+n*block] = 0.  # unused, nonexistent western face
                        continue
                    x, y = geo['rightix'][cell_partner], geo['rightiy'][cell_partner]
                    if not (-1 <= x <= nx and -1 <= y <= ny):
                        raise ValueError('Reflected physical face is outside geometry')
                    j = x+1+(nx+2)*(y+1)
                else:
                    j = cell_partner
                values[i+n*block] = .5*(old[i+n*block]+sign*old[j+n*block])
    for name in ['na', 'ne', 'te', 'ti']:
        if min(data[name]) <= 0:
            raise ValueError('Nonpositive '+name)
    after = totals(data, nv)
    inventory_errors = [abs(a-b)/abs(a) for a, b in zip(before[0], after[0])]
    ee_error, ie_error = abs(before[1]-after[1])/before[1], abs(before[2]-after[2])/before[2]
    if max(inventory_errors+[ee_error, ie_error]) > 1e-12:
        raise ValueError('Projection conservation check failed')
    write_state(destination, version, fields)
    report = {'proton_mass_kg': proton_mass, 'input': str(source), 'output': str(destination), 'physical_cells': len(physical),
              'species_inventory_relative_errors': inventory_errors,
              'electron_thermal_energy_relative_error': ee_error,
              'ion_thermal_plus_parallel_kinetic_energy_relative_error': ie_error,
              'parallel_momentum': 'Projected onto odd reflection parity; forbidden even part removed',
              'cache_face_fields_projected': face_fields,
              'minimum_na': min(data['na']), 'minimum_ne': min(data['ne']),
              'minimum_te_J': min(data['te']), 'minimum_ti_J': min(data['ti']),
              'neutral_restart': 'Not carried; cold neutral equilibration required'}
    destination.with_name('state_projection.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--old-b2', type=Path, required=True)
    parser.add_argument('--new-b2', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--proton-mass', type=float, default=1.67262192369e-27,
                        help='Match the linked solver constant; EFGW 0.1.1 uses CODATA 2018')
    args = parser.parse_args()
    project(args.source, args.old_b2, args.new_b2, args.output, args.proton_mass)
