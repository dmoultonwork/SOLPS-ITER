#!/usr/bin/env python3
"""Find the mirror pairs of the additional surfaces (block 3b) of a formatted Eirene input file (fort.1).

usage: make_surface_pairing.py [--tolerance CM] fort.1 pairing.txt

Two surfaces are mirror partners when the end points of one are the mirror images (Z -> -Z, same R) of the end points of the
other, in reversed order (mirroring a line reverses its direction, so the two surfaces then face the same way), to within
--tolerance (default 0.01 cm, Eirene lengths are in cm). A surface whose own end points are mirror images of each other
(typically a surface crossing the midplane) is paired with itself. The result is written in the format that
symmetrise_eirene_inputs.py reads. The script stops without writing a file if some surface has no partner, or more than
one, and lists them: such surfaces need a decision by hand. It reads only block 3b and does not compare the surface
physics; symmetrise_eirene_inputs.py does that when it uses the pairing.
"""
import argparse
import math
from pathlib import Path


def surfaces(path):
    lines = Path(path).read_text().split('\n')
    start = next(k for k, l in enumerate(lines) if l.startswith('*** 3b.'))
    count = int(lines[start+1])
    k, out = start+2, []
    for s in range(1, count+1):
        if not lines[k].startswith('*'):
            raise ValueError(f'Surface {s}: unexpected line {k+1}')
        c = [lines[k+3][12*m:12*m+12] for m in range(6)]
        out.append(((float(c[0]), float(c[1])), (float(c[3]), float(c[4]))))
        k += 4
        if not lines[k].startswith('*'):
            k += 1
    return out


def find_pairs(surf, tol):
    def near(p, q):
        return math.dist(p, q) <= tol
    pairs, problems = {}, []
    for i, (p1, p2) in enumerate(surf, 1):
        m1, m2 = (p1[0], -p1[1]), (p2[0], -p2[1])
        candidates = [j for j, (q1, q2) in enumerate(surf, 1) if near(q1, m2) and near(q2, m1)]
        if len(candidates) == 1:
            pairs[i] = candidates[0]
        else:
            reason = 'no partner' if not candidates else f'several partners {candidates}'
            same = [j for j, (q1, q2) in enumerate(surf, 1) if near(q1, m1) and near(q2, m2)]
            if same and not candidates:
                reason += f' (surface(s) {same} match with the end points in the same order, not reversed)'
            problems.append((i, reason))
    for i, j in pairs.items():
        if pairs.get(j) != i:
            problems.append((i, f'partner {j} does not point back to it'))
    return pairs, problems


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fort1', type=Path)
    p.add_argument('pairing', type=Path, help='output file (must not exist)')
    p.add_argument('--tolerance', type=float, default=0.01, help='end point tolerance in cm (default 0.01)')
    a = p.parse_args()
    if a.pairing.exists():
        raise SystemExit(f'{a.pairing} exists; not overwriting')
    surf = surfaces(a.fort1)
    pairs, problems = find_pairs(surf, a.tolerance)
    if problems:
        print(f'{len(surf)} additional surfaces; {len(problems)} without a unique mirror partner (no file written):')
        for i, reason in sorted(set(problems)):
            print(f'  surface {i}: ends {surf[i-1]}: {reason}')
        raise SystemExit(1)
    own = sum(1 for i, j in pairs.items() if i == j)
    lines = ['# mirror partner of each additional surface, from make_surface_pairing.py (tolerance %g cm)' % a.tolerance]
    lines += [f'{i} {pairs[i]}' for i in sorted(pairs)]
    a.pairing.write_text('\n'.join(lines) + '\n')
    print(f'{len(surf)} additional surfaces: {(len(surf)-own)//2} mirror pairs, {own} self-mirror; wrote {a.pairing}')
