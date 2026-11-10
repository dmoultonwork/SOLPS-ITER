#!/usr/bin/env python3
"""Check that a B2.5 grid is a connected double null and find its up-down mirror blocks (read-only).

usage: find_reflection_blocks.py b2fgmtry [--tolerance 0.25]

Reads the topology arrays of the B2.5 geometry file. Poloidal blocks start at cells without a left neighbour
(leftix = -2) and end at cells without a right neighbour (rightix = nx+1); every block includes its two guard
cells, so the blocks are given as inclusive LO:HI pairs in the form used by symmetrise_carre.py and
audit_b2_reflection.py, together covering -1..nx. The topology is classified from the number of cuts and the
radial extent of the cuts (B2.5 conventions); for an up-down symmetric connected double null the grid has two
blocks, each mapped onto itself by x -> LO+HI-x with Z -> -Z. The cell centres of each block are then compared
with their mirror images: the largest distance is reported in metres and in units of the local cell size.
The classification rules follow the B2.5 cut conventions and were checked on STEP connected double null grids and on grids
without cuts; no single-null or disconnected double-null grid was available to test them.
The exit status is 0 only for a connected double null whose blocks are mirror images within the tolerance
(default 0.25 cell sizes; a grid that is nearly symmetric is the case the symmetrisation tools are meant for).
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_mesh_symmetry import cf_fields


def analyse(path, tolerance):
    f = cf_fields(Path(path), {'nx,ny', 'nncut', 'leftcut', 'rightcut', 'topcut', 'bottomcut', 'leftix', 'rightix',
                               'crx', 'cry', 'hx', 'hy'})
    nx, ny = f['nx,ny']
    W, H = nx + 2, ny + 2
    n = W * H
    idx = lambda ix, iy: (ix + 1) + W * (iy + 1)
    ncut = f['nncut'][0]
    cuts = {k: f[k][:ncut] for k in ('leftcut', 'rightcut', 'topcut', 'bottomcut')}
    print(f'grid nx={nx} ny={ny}; cuts: {ncut}')
    for i in range(ncut):
        print(f'  cut {i+1}: x {cuts["leftcut"][i]}..{cuts["rightcut"][i]}, rows {cuts["bottomcut"][i]+1}..{cuts["topcut"][i]}')
    # poloidal blocks, taken in the outermost radial row and checked in all rows
    def blocks_in_row(iy):
        starts = [ix for ix in range(-1, nx + 1) if f['leftix'][idx(ix, iy)] == -2]
        ends = [ix for ix in range(-1, nx + 1) if f['rightix'][idx(ix, iy)] == nx + 1]
        return starts, ends
    starts, ends = blocks_in_row(ny - 1)
    consistent = all(blocks_in_row(iy) == (starts, ends) for iy in range(-1, ny + 1))
    blocks = list(zip(starts, ends))
    ok = len(starts) == len(ends) and all(a < b for a, b in blocks) and blocks[0][0] == -1 and blocks[-1][1] == nx \
        and all(blocks[k][1] + 1 == blocks[k + 1][0] for k in range(len(blocks) - 1))
    print('poloidal blocks (inclusive guard-cell endpoints):', ', '.join(f'{a}:{b}' for a, b in blocks) if blocks else 'none',
          '' if consistent else '(NOT the same in every radial row)')
    verdict = []
    if not ok or not consistent:
        verdict.append('block structure not recognised')
    # topology
    if ncut == 1 and len(blocks) == 1:
        topo = 'single null'
    elif ncut == 2 and len(blocks) == 2 and cuts['topcut'][0] == cuts['topcut'][1] and cuts['bottomcut'][0] == cuts['bottomcut'][1]:
        topo = 'connected double null (both X-points on the same flux surface)'
    elif ncut == 2 and len(blocks) == 2:
        topo = 'disconnected double null (X-points on different flux surfaces)'
    elif ncut == 0:
        topo = 'no X-point cuts (not a divertor grid)'
    else:
        topo = 'unrecognised topology'
    print('topology:', topo)
    if not topo.startswith('connected double null'):
        verdict.append('not a connected double null')
    # mirror test of the cell centres of each block (physical cells only)
    def centre(ix, iy):
        i = idx(ix, iy)
        return sum(f['crx'][i + n * c] for c in range(4)) / 4, sum(f['cry'][i + n * c] for c in range(4)) / 4
    worst = 0.0
    for a, b in blocks:
        wb = 0.0
        for iy in range(0, ny):
            for ix in range(a + 1, b):
                r1, z1 = centre(ix, iy)
                r2, z2 = centre(a + b - ix, iy)
                size = min(f['hx'][idx(ix, iy)], f['hy'][idx(ix, iy)])
                d = ((r1 - r2) ** 2 + (z1 + z2) ** 2) ** 0.5
                wb = max(wb, d / size if size > 0 else 0.0)
                if d / size > worst:
                    worst, wd = d / size, d
        print(f'  block {a}:{b}: largest mismatch of mirrored cell centres {wb:.3g} cell sizes')
    print(f'largest mismatch overall: {worst:.3g} cell sizes (tolerance {tolerance})')
    if worst > tolerance:
        verdict.append('not up-down mirror symmetric within the tolerance')
    if blocks:
        spec = ' '.join(f'{a}:{b}' for a, b in blocks)
        print('suggested:  audit_b2_reflection.py', path, spec)
        print('            symmetrise_carre.py ' + ' '.join(f'--block {a}:{b}' for a, b in blocks) + ' source.geo destination.geo')
    print('RESULT:', 'OK' if not verdict else 'not suitable: ' + '; '.join(verdict))
    return not verdict


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('b2fgmtry', type=Path)
    p.add_argument('--tolerance', type=float, default=0.25, help='largest allowed centre mismatch in cell sizes')
    a = p.parse_args()
    sys.exit(0 if analyse(a.b2fgmtry, a.tolerance) else 1)
