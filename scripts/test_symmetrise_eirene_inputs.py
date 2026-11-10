"""Synthetic fixtures for symmetrise_eirene_inputs.py (formatted fort.1 with two mirror surfaces and one surface source)."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('sei', Path(__file__).with_name('symmetrise_eirene_inputs.py'))
sei = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sei)


def fort1(upper_xy, upper_source, label_upper=-3):
    f = lambda *v: ''.join('%12.5E' % x for x in v)
    return '\n'.join([
        '*** 3b. Data for additional surfaces',
        '     2',
        '*    1 :    1', ' 2.00000E+00 1.00000E+00 1.00000E+00 1.00000E-05',
        '     0     0     0     0     0     1     0     0     0    -4',
        f(244.0, -642.291, -1e20, 243.0, -641.593, 1e20),
        '*    2 :    2', ' 2.00000E+00 1.00000E+00 1.00000E+00 1.00000E-05',
        '     0     0     0     0     0     1     0     0     0    %2d' % label_upper,
        f(*upper_xy),
        '*** 4. Data for species and atomic physics module',
        '*** 7. Data for primary sources, nstrai strata',
        '     1',
        '*  32 : Gas puffing source  2 : Ar', 'FFFFF', '  9999    -1   116     1     0',
        ' 1.00000E+00 0.00000E+00     0     0     0     0     0     0     0     0',
        'TFFFF', '     2', 'FFTFF', '     2',
        '    41     0     1     0     0     0     0     0     0',
        ' 3.00000E+00 2.20000E+02 0.00000E+00 0.00000E+00 1.00000E+03',
        '     0     0     0     0     0     0     0',
        f(243.0, 244.0, -641.593, -642.291, 0.0, 360.0),
        '    42     0     2     0     0     0     0     0     0',
        ' 3.00000E+00 2.20000E+02 0.00000E+00 0.00000E+00 1.00000E+03',
        '     0     0     0     0     0     0     0',
        f(*upper_source),
        '*** 8. Additional data for specific zones', ''])


class TestSymmetriseEireneInputs(unittest.TestCase):
    def run_case(self, text, **kw):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d/'in').write_text(text)
            (d/'pairs').write_text('# pairs\n1 2\n2 1\n')
            sei.prepare(d/'in', d/'pairs', d/'out', kw.get('swaps', [(-4, -3)]))
            return (d/'out').read_text()

    def test_mirrors_surface_and_source_exactly(self):
        out = self.run_case(fort1((243.0, 641.593, -1e20, 243.999, 642.290) + (1e20,),
                                  (243.0, 243.999, 641.593, 642.290, 0.0, 360.0)))
        lines = out.split('\n')
        self.assertEqual(lines[9], ''.join('%12.5E' % x for x in (243.0, 641.593, -1e20, 244.0, 642.291, 1e20)))
        self.assertIn(''.join('%12.5E' % x for x in (243.0, 244.0, 641.593, 642.291, 0.0, 360.0)), lines)
        self.assertEqual(self.run_case(out), out)          # already symmetric: nothing changes

    def test_stops_on_physics_mismatch(self):
        with self.assertRaises(ValueError):
            self.run_case(fort1((243.0, 641.593, -1e20, 244.0, 642.291, 1e20), (243.0, 244.0, 641.593, 642.291, 0.0, 360.0), label_upper=-9))


if __name__ == '__main__':
    unittest.main()
