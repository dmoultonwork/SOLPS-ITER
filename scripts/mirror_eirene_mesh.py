#!/usr/bin/env python3
"""Build a symmetric candidate Eirene mesh from a verified lower half.

Moves plasma vertices with the B2 quadrilateral map, clips vacuum triangles at
Z=0, reflects the half mesh, and reconstructs reciprocal edge neighbours. This
does not certify surface/source physics: those inputs must be audited separately.
"""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

from audit_mesh_symmetry import cf_fields


def cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def sub(a, b):
    return (a[0]-b[0], a[1]-b[1])


def area(points):
    a, b, c = points
    return .5*cross(sub(b, a), sub(c, a))


def toroidal_volume(points):
    return 2*math.pi*abs(area(points))*sum(p[0] for p in points)/3/1e6


def interpolate_quad(point, old, new):
    # Exact corner/edge treatment keeps the B2 coupling boundary conforming.
    tolerance = 2e-6  # cm: rounding in the original fort.33
    for a, b in zip(old, new):
        if math.dist(point, a) <= tolerance:
            return b
    for i, j in [(0, 1), (1, 3), (3, 2), (2, 0)]:
        delta = sub(old[j], old[i])
        length2 = sum(v*v for v in delta)
        if not length2:
            continue
        t = sum(a*b for a, b in zip(sub(point, old[i]), delta))/length2
        projected = tuple(old[i][k]+t*delta[k] for k in range(2))
        if -.0001 <= t <= 1.0001 and math.dist(projected, point) <= tolerance:
            t = min(1., max(0., t))
            return tuple((1-t)*new[i][k]+t*new[j][k] for k in range(2))
    u = v = .5
    for _ in range(30):
        weights = [(1-u)*(1-v), u*(1-v), (1-u)*v, u*v]
        mapped = tuple(sum(w*p[k] for w, p in zip(weights, old)) for k in range(2))
        du = tuple((1-v)*(old[1][k]-old[0][k])+v*(old[3][k]-old[2][k]) for k in range(2))
        dv = tuple((1-u)*(old[2][k]-old[0][k])+u*(old[3][k]-old[1][k]) for k in range(2))
        residual = sub(mapped, point)
        determinant = cross(du, dv)
        if abs(determinant) < 1e-25:
            raise ValueError('Degenerate B2 quadrilateral')
        step_u, step_v = cross(residual, dv)/determinant, cross(du, residual)/determinant
        u -= step_u
        v -= step_v
        if max(abs(step_u), abs(step_v)) < 1e-12:
            break
    if not (-.001 <= u <= 1.001 and -.001 <= v <= 1.001):
        raise ValueError(f'Eirene vertex outside assigned B2 quadrilateral: {u}, {v}')
    weights = [(1-u)*(1-v), u*(1-v), (1-u)*v, u*v]
    return tuple(sum(w*p[k] for w, p in zip(weights, new)) for k in range(2))


def neighbours(triangles):
    edges = defaultdict(list)
    for i, triangle in enumerate(triangles):
        for side in range(3):
            edges[tuple(sorted((triangle[side], triangle[(side+1) % 3])))].append((i, side))
    result = [[[0, 0, -1] for _ in range(3)] for t in triangles]
    boundary = 0
    for edge, entries in edges.items():
        if len(entries) == 1:
            boundary += 1
        elif len(entries) == 2:
            (i, side), (j, other) = entries
            if i == j:
                raise ValueError('Degenerate triangle edge')
            result[i][side] = [j+1, other+1, 0]
            result[j][other] = [i+1, side+1, 0]
        else:
            raise ValueError(f'Non-manifold edge {edge}: {entries}')
    return result, boundary


def b2_x(ix):
    # Verified extended poloidal indexing from original fort.30/35.
    for lo, hi, offset in [(1, 12, 1), (14, 61, 2), (63, 74, 3),
                           (78, 113, 4), (115, 162, 5), (164, 199, 6)]:
        if lo <= ix <= hi:
            return ix-offset
    raise ValueError(f'Unexpected plasma coarse index {ix}')


def reflected_coarse(ix, iy):
    if ix <= 0 and iy <= 0:
        return ix, iy
    if 1 <= ix <= 74:
        return 75-ix, iy
    if 78 <= ix <= 199:
        return 277-ix, iy
    raise ValueError(f'Cannot reflect coarse cell {ix}, {iy}')


def generate(original, old_geometry, new_geometry, output):
    f = cf_fields(old_geometry, {'nx,ny', 'crx', 'cry', 'vol'})
    g = cf_fields(new_geometry, {'nx,ny', 'crx', 'cry', 'vol'})
    if f['nx,ny'] != [194, 32] or g['nx,ny'] != f['nx,ny']:
        raise ValueError('This verified coarse mapping applies only to STEP 194x32')
    nx, ny = f['nx,ny']
    count = (nx+2)*(ny+2)

    def corners(data, ix, iy):
        index = b2_x(ix)+1+(nx+2)*iy  # IYTRI is one-based.
        return [(100*data['crx'][index+k*count], 100*data['cry'][index+k*count]) for k in range(4)]

    tokens = (original/'fort.33').read_text().split()
    n = int(tokens[0])
    if len(tokens) != 1+2*n:
        raise ValueError('Bad original vertex layout')
    points = list(zip(map(float, tokens[1:n+1]), map(float, tokens[n+1:])))
    triangles = [tuple(v-1 for v in map(int, line.split()[1:]))
                 for line in (original/'fort.34').read_text().splitlines()[1:]]
    records = [list(map(int, line.split())) for line in (original/'fort.35').read_text().splitlines()[1:]]
    if len(triangles) != len(records):
        raise ValueError('Triangle/record count mismatch')
    original_neighbours, old_boundary = neighbours(triangles)
    for i, row in enumerate(records):
        if row[0] != i+1:
            raise ValueError('Original triangle records are not ordered')
        for side in range(3):
            if row[1+3*side:3+3*side] != original_neighbours[i][side][:2]:
                raise ValueError(f'Original edge topology mismatch at triangle {i+1}')
    proposals = defaultdict(list)
    old_cells, new_cells = {}, {}
    for triangle, row in zip(triangles, records):
        ix, iy = row[-2:]
        if ix <= 0 or iy <= 0:
            continue
        key = ix, iy
        if key not in old_cells:
            old_cells[key], new_cells[key] = corners(f, ix, iy), corners(g, ix, iy)
        for v in triangle:
            proposals[v].append(interpolate_quad(points[v], old_cells[key], new_cells[key]))
    moved = points[:]
    conflict = 0.
    for v, options in proposals.items():
        moved[v] = options[0]
        conflict = max(conflict, max(math.dist(options[0], q) for q in options))
    if conflict > 1e-7:
        raise ValueError(f'Inconsistent shared B2 vertex mapping: {conflict} cm')
    max_shift = max(math.dist(p, q) for p, q in zip(points, moved))
    # The raw file's midplane is displaced by less than 0.001 cm. Explicitly
    # snap those shared vertices, then clip the genuinely crossing vacuum cells.
    moved = [(r, 0. if abs(z) < .001 else z) for r, z in moved]
    intersections = {}

    def intersection(a, b):
        if moved[a][1] == 0.:
            return a
        if moved[b][1] == 0.:
            return b
        key = tuple(sorted((a, b)))
        if key not in intersections:
            pa, pb = moved[a], moved[b]
            t = -pa[1]/(pb[1]-pa[1])
            intersections[key] = len(moved)
            moved.append((pa[0]+t*(pb[0]-pa[0]), 0.))
        return intersections[key]

    half, metadata, origin = [], [], []
    clipped = 0
    for i, (triangle, row) in enumerate(zip(triangles, records)):
        polygon = []
        for a, b in zip(triangle, triangle[1:]+triangle[:1]):
            inside_a, inside_b = moved[a][1] <= 0, moved[b][1] <= 0
            if inside_a:
                polygon.append(a)
            if inside_a != inside_b:
                polygon.append(intersection(a, b))
        polygon = list(dict.fromkeys(polygon))
        if len(polygon) < 3:
            continue
        if len(polygon) > 3:
            clipped += 1
        for k in range(1, len(polygon)-1):
            t = (polygon[0], polygon[k], polygon[k+1])
            signed_area = area([moved[v] for v in t])
            if abs(signed_area) < 1e-18:
                continue
            if signed_area < 0:
                t = (t[0], t[2], t[1])
            half.append(t)
            metadata.append(tuple(row[-2:]))
            origin.append(i)
    used = sorted({v for t in half for v in t})
    remap = {v: i for i, v in enumerate(used)}
    full_points = [moved[v] for v in used]
    vertex_map = []
    for i, (r, z) in enumerate(full_points[:]):
        if z == 0.:
            vertex_map.append(i)
        else:
            vertex_map.append(len(full_points))
            full_points.append((r, -z))
    vertex_map.extend(i for i in range(len(used)) if full_points[i][1] != 0.)
    if any(vertex_map[vertex_map[i]] != i for i in range(len(full_points))):
        raise ValueError('Vertex reflection is not involutive')
    half = [tuple(remap[v] for v in t) for t in half]
    full_triangles = half+[tuple(vertex_map[v] for v in (t[0], t[2], t[1])) for t in half]
    full_metadata = metadata+[reflected_coarse(*cell) for cell in metadata]
    new_neighbours, boundary = neighbours(full_triangles)
    volumes = defaultdict(float)
    for t, cell in zip(full_triangles, full_metadata):
        p = [full_points[v] for v in t]
        if area(p) <= 0:
            raise ValueError('Nonpositive new triangle area')
        if cell[0] > 0 and cell[1] > 0:
            volumes[cell] += toroidal_volume(p)
    # The mirrored mesh must cover the same B2 cells as the original, and no cell may be worse than the original Eirene
    # mesh already was against the original B2 volume (a few cells of a given mesh can differ from the B2 volume by
    # ~1e-5 for reasons unrelated to symmetrisation; the mirroring preserves that and must not add to it).
    original_volumes = defaultdict(float)
    for triangle, row in zip(triangles, records):
        if row[-2] > 0 and row[-1] > 0:
            original_volumes[tuple(row[-2:])] += toroidal_volume([points[v] for v in triangle])
    if set(volumes) != set(original_volumes):
        raise ValueError(f'B2 coupling coverage mismatch: {len(volumes)} cells now, {len(original_volumes)} originally')
    mismatch = inherited = 0.
    worse = []
    for (ix, iy), value in volumes.items():
        index = b2_x(ix)+1+(nx+2)*iy
        err = abs(value-g['vol'][index])/g['vol'][index]
        old_err = abs(original_volumes[(ix, iy)]-f['vol'][index])/f['vol'][index]
        mismatch, inherited = max(mismatch, err), max(inherited, old_err)
        if err > old_err+1e-8:
            worse.append((err, old_err, ix, iy))
    if worse:
        worse.sort(reverse=True)
        raise ValueError(f'B2 coupling volume mismatch: {len(worse)} cells are worse than the original mesh; worst '
                         f'(new error, original error, ix, iy) = {worse[0]}')
    output.mkdir(parents=True, exist_ok=True)
    with (output/'fort.33').open('w') as stream:
        stream.write(f'{len(full_points):12d}\n')
        for axis in range(2):
            for i in range(0, len(full_points), 4):
                stream.write(' '.join(f'{p[axis]:.17e}' for p in full_points[i:i+4])+'\n')
    with (output/'fort.34').open('w') as stream:
        stream.write(f'{len(full_triangles):12d}\n')
        for i, t in enumerate(full_triangles):
            stream.write(' '.join(map(str, [i+1]+[v+1 for v in t]))+'\n')
    with (output/'fort.35').open('w') as stream:
        stream.write(f'{len(full_triangles):12d}\n')
        for i, (edges, cell) in enumerate(zip(new_neighbours, full_metadata)):
            stream.write(' '.join(map(str, [i+1]+[v for edge in edges for v in edge]+list(cell)))+'\n')
    report = {'scope': 'Candidate geometry/connectivity; surface/source physics still require verification',
              'original_vertices': n, 'original_triangles': len(triangles),
              'vertices': len(full_points), 'triangles': len(full_triangles),
              'old_boundary_edges': old_boundary, 'new_boundary_edges': boundary,
              'midplane_intersections': len(intersections), 'split_lower_quadrilaterals': clipped,
              'maximum_original_vertex_shift_cm': max_shift,
              'maximum_shared_vertex_disagreement_cm': conflict,
              'coupled_b2_cells': len(volumes), 'maximum_coupling_volume_relative_error': mismatch,
              'maximum_original_mesh_volume_relative_error': inherited,
              'vertex_reflection_zero_based': vertex_map,
              'triangle_reflection_zero_based': list(range(len(half), 2*len(half)))+list(range(len(half)))}
    (output/'mesh_generation.json').write_text(json.dumps(report, indent=2)+'\n')
    # Reflection map of the Eirene cells (zero-based, one index per line), read by the neutral-restart tools.
    # Index (zero-based) of the original triangle each new triangle comes from: the lower half keeps its triangles (a
    # clipped triangle gives several), the upper half is the mirror image of the lower half and takes the same origin.
    (output/'triangle_origin.txt').write_text(''.join('%d\n' % i for i in origin+origin))
    (output/'triangle_reflection.txt').write_text(''.join('%d\n' % i for i in report['triangle_reflection_zero_based']))
    print(json.dumps({k: v for k, v in report.items() if not isinstance(v, list)}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--old-b2', type=Path, required=True)
    parser.add_argument('--new-b2', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    generate(args.original, args.old_b2, args.new_b2, args.output)
