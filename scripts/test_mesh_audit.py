"""Small analytic fixtures for reflection, ghost exclusion and mesh connectivity."""
import importlib.util
from pathlib import Path
import random
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('audit', Path(__file__).with_name('audit_mesh_symmetry.py'))
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class MeshAuditTests(unittest.TestCase):
    def test_nearest_matches_brute_force(self):
        rng = random.Random(19)
        points = [(rng.uniform(-5, 5), rng.uniform(-5, 5)) for _ in range(301)]
        kd = audit.tree(points)
        for _ in range(100):
            query = (rng.uniform(-6, 6), rng.uniform(-6, 6))
            expected = min((sum((a-b)**2 for a, b in zip(p, query)), i)
                           for i, p in enumerate(points))
            self.assertEqual(audit.nearest(points, kd, query), expected)

    def test_b2_symmetric_rectangle_and_volume_error(self):
        coordinates = [(x, y) for y in range(-1, 3) for x in range(-1, 3)]
        fields = {'nx,ny': ('int', [2, 2])}
        fields['crx'] = ('real', [2+y+(k//2) for k in range(4) for x, y in coordinates])
        fields['cry'] = ('real', [x-1+(k % 2) for k in range(4) for x, y in coordinates])
        for name in ['vol', 'hx', 'hy']:
            fields[name] = ('real', [1.]*16)
        for name, dx, dy in [('left', -1, 0), ('right', 1, 0), ('top', 0, 1), ('bottom', 0, -1)]:
            fields[name+'ix'] = ('int', [x+dx for x, y in coordinates])
            fields[name+'iy'] = ('int', [y+dy for x, y in coordinates])

        def write(path):
            path.write_text('VERSION\n'+''.join(f'*cf: {kind} {len(values)} {name}\n'+
                            ' '.join(map(str, values))+'\n'
                            for name, (kind, values) in fields.items()))

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'b2fgmtry'
            write(path)
            result = audit.b2_audit(path, 0.)['logical_reflection']
            self.assertEqual(result['physical_cell_count'], 4)
            self.assertEqual(result['connectivity_failures_by_direction'], {})
            self.assertEqual(result['corner_mismatch_native_units']['maximum'], 0.)
            self.assertEqual(result['metric_relative_mismatch']['vol']['maximum'], 0.)
            fields['vol'][1][5] = 2.
            write(path)
            result = audit.b2_audit(path, 0.)['logical_reflection']
            self.assertEqual(result['metric_relative_mismatch']['vol']['maximum'], .5)

    def test_eirene_symmetric_fan_and_missing_reflected_triangle(self):
        with tempfile.TemporaryDirectory() as tmp:
            vertices, triangles = Path(tmp)/'fort.33', Path(tmp)/'fort.34'
            vertices.write_text('5\n1 2 1 2 1.5\n-1 -1 1 1 0\n')
            triangles.write_text('4\n1 1 2 5\n2 2 4 5\n3 4 3 5\n4 3 1 5\n')
            result = audit.eirene_audit(vertices, triangles, 0.)
            self.assertEqual(result['distance_native_units']['maximum'], 0.)
            self.assertEqual(result['triangles_without_candidate_reflected_connectivity'], 0)
            triangles.write_text('3\n1 1 2 5\n2 2 4 5\n3 3 1 5\n')
            result = audit.eirene_audit(vertices, triangles, 0.)
            self.assertEqual(result['triangles_without_candidate_reflected_connectivity'], 1)


if __name__ == '__main__':
    unittest.main()
