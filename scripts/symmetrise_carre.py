#!/usr/bin/env python3
"""Create a reflected CARRE .geo input; regenerate B2 metrics with native b2ag.

This does not edit a b2fgmtry or assert that Eirene or the boundary data are symmetric.
The logical B2 reflection must have been checked before supplying the block endpoints.
"""
import argparse
import json
from pathlib import Path


def symmetrise(source, destination, blocks):
    lines = source.read_text().splitlines()
    nx, ny = map(int, lines[0].split())
    count = (nx+2)*(ny+2)
    if len(lines) != 1+2*count:
        raise ValueError('Expected CARRE coordinates/fields plus flux block')
    rows, flux = [], []
    for start, target, width in [(1, rows, 12), (1+count, flux, 8)]:
        for i, line in enumerate(lines[start:start+count]):
            tokens = line.split()
            x, y = map(int, tokens[:2])
            if x+(nx+2)*y != i or len(tokens) != 2+width:
                raise ValueError('Unexpected row layout')
            target.append(list(map(float, tokens[2:])))
    xmap = {x+1: lo+hi-x+1 for lo, hi in blocks for x in range(lo, hi+1)}
    if set(xmap) != set(range(nx+2)) or any(xmap[xmap[x]] != x for x in xmap):
        raise ValueError('Reflection blocks must cover every column and be involutive')
    # Centre followed by four corners; reverse the logical poloidal direction.
    point_partner = [0, 2, 1, 4, 3]
    flux_partner = [1, 0, 3, 2, 5, 4, 7, 6]
    max_change = 0.
    for y in range(ny+2):
        for x in range(nx+2):
            partner = xmap[x]
            if x >= partner:
                continue
            i, j = x+(nx+2)*y, partner+(nx+2)*y
            old_a, old_b = rows[i][:], rows[j][:]
            for point, other in enumerate(point_partner):
                r = .5*(old_a[2*point]+old_b[2*other])
                z = .5*(old_a[2*point+1]-old_b[2*other+1])
                rows[i][2*point:2*point+2] = [r, z]
                rows[j][2*other:2*other+2] = [r, -z]
            # Bp and Bt in CARRE are signed components in the local basis.
            for k in [10, 11]:
                rows[i][k] = rows[j][k] = .5*(old_a[k]+old_b[k])
            for k in range(10):
                max_change = max(max_change, abs(rows[i][k]-old_a[k]), abs(rows[j][k]-old_b[k]))
            old_a, old_b = flux[i][:], flux[j][:]
            for k, other in enumerate(flux_partner):
                flux[i][k] = flux[j][other] = .5*(old_a[k]+old_b[other])
    # A shared midplane corner can be reflected onto itself through an adjacent
    # cell. Averaging equal coordinates yields exactly zero for its height.
    text = f'{nx:12d}{ny:12d}\n'
    for data in [rows, flux]:
        text += ''.join(f'{i%(nx+2):5d}{i//(nx+2):5d} '+
                        ' '.join(f'{v:.17e}' for v in row)+'\n' for i, row in enumerate(data))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text)
    return {'input': str(source), 'output': str(destination), 'shape': [nx, ny],
            'blocks_including_guards': blocks, 'maximum_coordinate_change_m': max_change,
            'metrics_regenerated': False, 'eirene_regenerated': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--block', action='append', required=True, help='Inclusive B2 guard endpoints, e.g. -1:72')
    args = parser.parse_args()
    blocks = [tuple(map(int, value.split(':'))) for value in args.block]
    print(json.dumps(symmetrise(args.source, args.destination, blocks), indent=2))
