#!/usr/bin/env python3
"""Mirror the additional-surface and surface-source coordinates of a formatted Eirene input file (fort.1).

usage: symmetrise_eirene_inputs.py [--swap A:B ...] fort.1 pairing.txt new_fort.1

pairing.txt lists the mirror partner of every additional surface (block 3b, numbered 1..N): one "i j" pair per line
(blank lines and lines starting with # are ignored); "i i" marks a surface that is its own mirror image. The mapping
must be a complete involution.

For every pair the surface in the upper half (Z > 0) gets exactly the mirror image of the end points of its partner in the
lower half (ends swapped, Z negated), and the same is done for the SORAD1..4 coordinates of the surface sources (strata with
a surface source, block 7). Coordinates are copied as text, so the mirror is exact in the file. All other parameters of
paired surfaces and sources must agree (the surface labels in the last integer of the surface line may differ by the
label swaps given with --swap, default -4:-3 and 5:6, which are those of the STEP case); otherwise the script stops.
Nothing else in the file is touched. It does not check the physics of the labels, the mesh or the plasma background:
those need separate checks. Neutral restart data are not handled here.
"""
import argparse
from pathlib import Path


def fields(line, n):
    """Fixed-width Eirene floats (12 characters each)."""
    return [line[12*k:12*k+12] for k in range(n)]


def neg(tok):
    if float(tok) == 0.0:
        return tok
    return ('-' + tok[1:]) if tok[0] == ' ' else (' ' + tok[1:])


def put(line, toks):
    return ''.join(toks) + line[12*len(toks):]


def parse_pairs(path, count):
    pairs = {}
    for raw in Path(path).read_text().splitlines():
        text = raw.split('#')[0].split()
        if text:
            i, j = int(text[0]), int(text[1])
            pairs[i] = j
    if sorted(pairs) != list(range(1, count+1)) or any(pairs[pairs[i]] != i for i in pairs):
        raise ValueError('The pairing must list every additional surface once and be a complete involution')
    return pairs


def prepare(source, pairing, destination, swaps):
    lines = Path(source).read_text().split('\n')
    label_map = {}
    for a, b in swaps:
        label_map[a] = b
        label_map[b] = a
    # --- block 3b: additional surfaces
    start = next(k for k, l in enumerate(lines) if l.startswith('*** 3b.'))
    count = int(lines[start+1])
    pairs = parse_pairs(pairing, count)
    k = start + 2
    surf = {}
    for s in range(1, count+1):
        if not lines[k].startswith('*'):
            raise ValueError(f'Surface {s}: unexpected line {k+1}')
        entry = {'par': k+1, 'ints': k+2, 'xy': k+3}
        k += 4
        if not lines[k].startswith('*'):
            k += 1                      # surface-model text line
        surf[s] = entry
    changes = []
    for i, j in sorted(pairs.items()):
        a, b = surf[i], surf[j]
        pa = lines[a['par']].split(); pb = lines[b['par']].split()
        ia = [int(x) for x in lines[a['ints']].split()]; ib = [int(x) for x in lines[b['ints']].split()]
        if pa != pb or ia[:-1] != ib[:-1] or label_map.get(ia[-1], ia[-1]) != ib[-1]:
            raise ValueError(f'Surface physics mismatch between {i} and {j}')
        ca = fields(lines[a['xy']], 6); cb = fields(lines[b['xy']], 6)
        if i == j:
            if ca[0] != cb[3] or float(ca[1]) != -float(cb[4]):
                raise ValueError(f'Surface {i} is its own mirror image but is not symmetric')
            continue
        if i > j:
            continue
        lower, upper = ((i, j) if float(ca[1]) + float(ca[4]) < 0 else (j, i))
        lo, up = surf[lower], surf[upper]
        cl = fields(lines[lo['xy']], 6); cu = fields(lines[up['xy']], 6)
        new = list(cu)
        new[0], new[1] = cl[3], neg(cl[4])        # upper P1 = mirror of lower P2
        new[3], new[4] = cl[0], neg(cl[1])        # upper P2 = mirror of lower P1
        if new != cu:
            changes.append((upper, lines[up['xy']], put(lines[up['xy']], new)))
            lines[up['xy']] = put(lines[up['xy']], new)
    # --- block 7: surface sources
    b7 = next(k for k, l in enumerate(lines) if l.startswith('*** 7.'))
    end7 = next(k for k in range(b7+1, len(lines)) if lines[k].startswith('*** '))
    heads = [k for k in range(b7, end7) if lines[k].startswith('*') and not lines[k].startswith('***')
             and lines[k].split(':')[0].strip('* ').isdigit()]
    source_changes = []
    for h in heads:
        flags = lines[h+6].strip()
        if len(flags) != 5 or set(flags) - set('TF'):
            raise ValueError(f'Stratum at line {h+1}: unexpected layout')
        if flags[2] != 'T':
            continue                              # not a surface source
        n = int(lines[h+7])
        rows = {}
        for r in range(n):
            base = h + 8 + 4*r
            ints = lines[base].split()
            rows[int(ints[2])] = base
        if len(rows) != n:
            raise ValueError(f'Stratum at line {h+1}: duplicate source surface')
        for i, base_i in sorted(rows.items()):
            j = pairs[i]
            if j not in rows:
                raise ValueError(f'Stratum at line {h+1}: source on surface {i} has no source on its mirror surface {j}')
            if i > j:
                continue
            base_j = rows[j]
            xa = fields(lines[base_i+3], 6); xb = fields(lines[base_j+3], 6)
            if i == j:
                if xa[0] != xa[1] or float(xa[2]) != -float(xa[3]):
                    raise ValueError(f'Source on surface {i} is its own mirror image but is not symmetric')
                continue
            ta = lines[base_i].split(); tb = lines[base_j].split()
            same = (ta[1:2] + ta[3:] == tb[1:2] + tb[3:] and lines[base_i+1] == lines[base_j+1] and lines[base_i+2] == lines[base_j+2]
                    and xa[4:] == xb[4:])
            if not same:
                raise ValueError(f'Source parameter mismatch between surfaces {i} and {j}')
            lower, upper = ((base_i, base_j) if float(xa[2]) + float(xa[3]) < 0 else (base_j, base_i))
            xl = fields(lines[lower+3], 6); xu = fields(lines[upper+3], 6)
            new = [xl[0], xl[1], neg(xl[2]), neg(xl[3])] + xu[4:]
            if new != xu:
                source_changes.append((int(lines[h].split(':')[0].strip('* ')), lines[upper+3], put(lines[upper+3], new)))
                lines[upper+3] = put(lines[upper+3], new)
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    destination.write_text('\n'.join(lines))
    print(f'{len(changes)} surface(s) and {len(source_changes)} source line(s) changed; wrote {destination}')
    for s, old, new in changes:
        print(f'  surface {s}:\n    {old.strip()}\n -> {new.strip()}')
    for st, old, new in source_changes:
        print(f'  stratum {st} source:\n    {old.strip()}\n -> {new.strip()}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('source', type=Path, help='formatted Eirene input file (fort.1)')
    p.add_argument('pairing', type=Path, help='text file with one "i j" mirror pair of additional surfaces per line')
    p.add_argument('destination', type=Path, help='new fort.1 (must not exist)')
    p.add_argument('--swap', action='append', default=None, metavar='A:B',
                   help='surface labels (last integer of the surface line) exchanged by the mirror, e.g. -4:-3; repeatable')
    a = p.parse_args()
    swaps = [tuple(int(v) for v in s.split(':')) for s in (a.swap or ['-4:-3', '5:6'])]
    prepare(a.source, a.pairing, a.destination, swaps)
