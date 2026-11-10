#!/usr/bin/env python3
"""Read-only coordinate and candidate-connectivity audit; no third-party packages.

Distances use each file's native units. B2 geometry normally uses metres and
Eirene fort.33 normally uses centimetres: specify/verify units before comparison.
Nearest reflected points are candidates, never proof of correct physical pairing.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sys


def cf_fields(path, wanted):
    fields = {}
    header = re.compile(r'^\*cf:\s+(int|real|char)\s+(\d+)\s+(\S+)')
    with path.open() as stream:
        for line in stream:
            match = header.match(line)
            if not match:
                continue
            kind, count, name = match.groups()
            count = int(count)
            if name not in wanted or kind == 'char':
                continue
            values = []
            while len(values) < count:
                line = next(stream)
                if line.startswith('*cf:'):
                    raise ValueError(f'Truncated field {name}')
                values.extend(line.replace('D', 'E').replace('d', 'e').split())
            if len(values) != count:
                raise ValueError(f'Incorrect count for {name}')
            fields[name] = [int(v) if kind == 'int' else float(v) for v in values]
    if wanted - fields.keys():
        raise ValueError(f'Missing fields: {wanted - fields.keys()}')
    return fields


def tree(points, indices=None, depth=0):
    if indices is None:
        indices = list(range(len(points)))
    if not indices:
        return None
    axis = depth % 2
    indices.sort(key=lambda i: (points[i][axis], i))
    mid = len(indices) // 2
    return (indices[mid], axis, tree(points, indices[:mid], depth+1),
            tree(points, indices[mid+1:], depth+1))


def nearest(points, node, query, best=(math.inf, -1)):
    if node is None:
        return best
    index, axis, left, right = node
    distance = sum((a-b)**2 for a, b in zip(points[index], query))
    best = min(best, (distance, index))
    delta = query[axis] - points[index][axis]
    close, far = (left, right) if delta <= 0 else (right, left)
    best = nearest(points, close, query, best)
    if delta*delta <= best[0]:
        best = nearest(points, far, query, best)
    return best


def statistics(values):
    values = sorted(values)
    return {'count': len(values), 'rms': math.sqrt(sum(x*x for x in values)/len(values)),
            'median': values[len(values)//2], 'p99': values[min(len(values)-1, int(.99*len(values)))],
            'maximum': values[-1]}


def reflection(points, plane):
    kd = tree(points)
    results = [nearest(points, kd, (r, 2*plane-z)) for r, z in points]
    mapping = [i for _, i in results]
    scale = max(1., max(abs(x) for p in points for x in p))
    tolerance = 128*sys.float_info.epsilon*scale
    distances = [math.sqrt(d) for d, _ in results]
    return mapping, {
        'distance_native_units': statistics(distances),
        'roundoff_tolerance_native_units': tolerance,
        'points_above_roundoff_tolerance': sum(d > tolerance for d in distances),
        'unique_candidate_partners': len(set(mapping)),
        'non_involutive_candidates': sum(mapping[j] != i for i, j in enumerate(mapping)),
    }


def b2_audit(path, plane):
    neighbours = {axis+suffix for axis in ['left', 'right', 'top', 'bottom'] for suffix in ['ix', 'iy']}
    f = cf_fields(path, {'nx,ny', 'crx', 'cry', 'vol', 'hx', 'hy'} | neighbours)
    nx, ny = f['nx,ny']
    count = (nx+2)*(ny+2)
    if len(f['crx']) != 4*count or len(f['cry']) != 4*count:
        raise ValueError('Unsupported B2 corner layout')
    def offset(x, y):
        return x+1+(nx+2)*(y+1)

    # Detect poloidal blocks by the external-neighbour sentinels in B2's own
    # connectivity. Each block includes a guard column at both ends.
    starts = [x for x in range(-1, nx+1)
              if all(f['leftix'][offset(x, y)] == -2 for y in range(ny))]
    ends = [x for x in range(-1, nx+1)
            if all(f['rightix'][offset(x, y)] == nx+1 for y in range(ny))]
    blocks = list(zip(starts, ends))
    if (not blocks or len(starts) != len(ends) or blocks[0][0] != -1 or
            blocks[-1][1] != nx or any(hi <= lo+1 for lo, hi in blocks) or
            any(b[0] != a[1]+1 for a, b in zip(blocks, blocks[1:]))):
        raise ValueError('Cannot infer complete poloidal blocks safely from B2 connectivity')
    xmap = {x: lo+hi-x for lo, hi in blocks for x in range(lo, hi+1)}
    xmap.update({-2: nx+1, nx+1: -2})
    physical_x = [x for lo, hi in blocks for x in range(lo+1, hi)]
    coordinates = [(x, y) for y in range(ny) for x in physical_x]
    indices = [offset(x, y) for x, y in coordinates]
    cells = [[(f['crx'][i+k*count], f['cry'][i+k*count]) for k in range(4)] for i in indices]
    centres = [(sum(p[0] for p in c)/4, sum(p[1] for p in c)/4) for c in cells]
    mapping, result = reflection(centres, plane)
    corner_errors = []
    volume_errors = []
    for i, j in enumerate(mapping):
        reflected = [(r, 2*plane-z) for r, z in cells[i]]
        corner_errors.append(min(max(math.dist(a, b) for a, b in zip(reflected, perm))
                                 for perm in itertools.permutations(cells[j])))
        va, vb = f['vol'][indices[i]], f['vol'][indices[j]]
        volume_errors.append(abs(va-vb)/max(abs(va), abs(vb), sys.float_info.min))
    result.update({'physical_shape': [nx, ny], 'guard_cells_checked': False,
                   'candidate_corner_mismatch_native_units': statistics(corner_errors),
                   'candidate_volume_relative_mismatch': statistics(volume_errors),
                   'topology_and_metrics_verified': False})
    # A second, explicit logical reflection is checked independently of nearest
    # neighbours. Only connectivity agreement makes it a usable pairing.
    topology_failures = {}
    corner_errors = []
    metric_errors = {name: [] for name in ['vol', 'hx', 'hy']}
    opposite = {'left': 'right', 'right': 'left', 'top': 'top', 'bottom': 'bottom'}
    for x, y in coordinates:
        i, j = offset(x, y), offset(xmap[x], y)
        for side, other in opposite.items():
            neighbour = (f[side+'ix'][i], f[side+'iy'][i])
            reflected = (xmap.get(neighbour[0]), neighbour[1])
            target = (f[other+'ix'][j], f[other+'iy'][j])
            if reflected != target:
                topology_failures[side] = topology_failures.get(side, 0)+1
        original = [(f['crx'][i+k*count], 2*plane-f['cry'][i+k*count]) for k in range(4)]
        partner = [(f['crx'][j+k*count], f['cry'][j+k*count]) for k in [1, 0, 3, 2]]
        corner_errors.append(max(math.dist(a, b) for a, b in zip(original, partner)))
        for name in metric_errors:
            va, vb = f[name][i], f[name][j]
            metric_errors[name].append(abs(va-vb)/max(abs(va), abs(vb), sys.float_info.min))
    result['logical_reflection'] = {
        'poloidal_blocks_including_guards': blocks,
        'physical_cell_count': len(coordinates),
        'connectivity_failures_by_direction': topology_failures,
        'corner_mismatch_native_units': statistics(corner_errors),
        'metric_relative_mismatch': {k: statistics(v) for k, v in metric_errors.items()},
        'scope': 'Physical-cell neighbour topology and corner parity; further field/guard/coupling checks required.'}
    return result


def eirene_audit(vertices, triangles, plane):
    tokens = vertices.read_text().replace('D', 'E').split()
    count = int(tokens[0])
    if len(tokens) != 1+2*count:
        raise ValueError('Unsupported fort.33 layout: expected count, R array, Z array')
    points = list(zip(map(float, tokens[1:1+count]), map(float, tokens[1+count:])))
    mapping, result = reflection(points, plane)
    lines = triangles.read_text().splitlines()
    ntri = int(lines[0])
    connectivity = {}
    for line in lines[1:]:
        row = list(map(int, line.split()))
        if len(row) != 4 or not 1 <= row[0] <= ntri:
            raise ValueError('Unsupported fort.34 layout')
        if any(not 1 <= v <= count for v in row[1:]):
            raise ValueError('Invalid triangle vertex')
        if row[0] in connectivity:
            raise ValueError('Duplicate triangle ID')
        connectivity[row[0]] = tuple(v-1 for v in row[1:])
    if len(connectivity) != ntri:
        raise ValueError('Triangle count mismatch')
    keys = {tuple(sorted(t)) for t in connectivity.values()}
    unmatched = sum(tuple(sorted(mapping[v] for v in t)) not in keys for t in connectivity.values())
    result.update({'vertices': count, 'triangles': ntri,
                   'triangles_without_candidate_reflected_connectivity': unmatched,
                   'candidate_connectivity_only': True,
                   'materials_surfaces_and_coupling_verified': False})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--b2-geometry', type=Path, required=True)
    parser.add_argument('--eirene-vertices', type=Path, required=True)
    parser.add_argument('--eirene-triangles', type=Path, required=True)
    parser.add_argument('--b2-plane', type=float, default=0.)
    parser.add_argument('--eirene-plane', type=float, default=0.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    paths = [args.b2_geometry, args.eirene_vertices, args.eirene_triangles]
    report = {'scope': 'Original mesh screening; candidate pairings require topology validation before editing.',
              'inputs': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
              'reflection_planes_native_units': {'b2': args.b2_plane, 'eirene': args.eirene_plane},
              'b2': b2_audit(args.b2_geometry, args.b2_plane),
              'eirene': eirene_audit(args.eirene_vertices, args.eirene_triangles, args.eirene_plane)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
