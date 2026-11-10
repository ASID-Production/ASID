import pathlib
import shutil
import tempfile
import unittest

from test_distance_loading import EXTENSION, _LOAD_EXTENSION, _run_isolated


class BaderApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not EXTENSION.is_file():
            raise unittest.SkipTest(f"Extension not found: {EXTENSION}")

    def test_invalid_arguments_preserve_python_exceptions(self):
        script = _LOAD_EXTENSION + r"""
cell = [4., 4., 4., 90., 90., 90.]
symm = ['x,y,z']
atoms = [(1, .25, .25, .25)]
cases = [(), (cell, symm, atoms), (cell, symm, atoms, 6., 0),
         (None, symm, atoms, 6.), ([4.], symm, atoms, 6.),
         (cell, None, atoms, 6.), (cell, [None], atoms, 6.),
         (cell, ['x'], atoms, 6.), (cell, ['x,y,1/0'], atoms, 6.),
         (cell, ['x,y,\ud800'], atoms, 6.),
         (cell, ['x,y,999999999999'], atoms, 6.),
         (cell, ['x,y,z\x00junk'], atoms, 6.),
         (cell, symm, None, 6.), (cell, symm, [], 6.),
         (cell, symm, [(1,)], 6.), (cell, symm, [[1, 0, 0, 0]], 6.),
         (cell, symm, [(0, 0, 0, 0)], 6.),
         (cell, symm, [(119, 0, 0, 0)], 6.),
         (cell, symm, [(1 << 200, 0, 0, 0)], 6.),
         (cell, symm, [(1, 'x', 0, 0)], 6.),
         (cell, symm, [(1, float('nan'), 0, 0)], 6.),
         (cell, symm, atoms, float('nan')),
         (cell, symm, atoms, float('inf')), (cell, symm, atoms, 0.),
         ([0.,4.,4.,90.,90.,90.], symm, atoms, 6.),
         ([4.,4.,4.,10.,10.,170.], symm, atoms, 6.)]
for arguments in cases:
    try:
        cpplib.FindCP(*arguments)
    except (TypeError, ValueError, OverflowError, UnicodeError):
        pass
    else:
        raise AssertionError(repr(arguments))
print('validated', len(cases))
"""
        result = _run_isolated(script, extension=EXTENSION)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('validated', result.stdout)

    def test_density_loading_errors_are_exceptions(self):
        for content in (None, '', '1\n1 .1 1 0\n', '1\n1 0 1 1\n2\n',
                        '1\n1 .1 1 1\nnan\n', '-1\n'):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as tmp:
                root = pathlib.Path(tmp)
                target = root / EXTENSION.name
                shutil.copy2(EXTENSION, target)
                if content is not None:
                    (root / 'ED.txt').write_text(content)
                script = _LOAD_EXTENSION + """
try:
    cpplib.FindCP([4.,4.,4.,90.,90.,90.], ['x,y,z'], [(1,.25,.25,.25)], 6.)
except RuntimeError as error:
    assert any(state in str(error) for state in ('FILE_NOT_FOUND', 'BAD_SIZE', 'BAD_CONTENT')), str(error)
else:
    raise AssertionError('invalid density accepted')
"""
                result = _run_isolated(script, extension=target, cwd=root)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_findcp_valid_crystal(self):
        with tempfile.TemporaryDirectory() as tmp:
            script = _LOAD_EXTENSION + """
cpplib.FindCP([4.,4.,4.,90.,90.,90.], ['x,y,z'], [(1,.25,.25,.25)], 6.)
assert pathlib.Path('CPs.pdb').is_file()
assert pathlib.Path('paths.pdb').is_file()
print('FindCP OK')
"""
            result = _run_isolated(script, extension=EXTENSION, cwd=tmp, timeout=60)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('FindCP OK', result.stdout)

    def test_voronoi_conversion_allocation_failures(self):
        script = _LOAD_EXTENSION + """
import _testcapi
arguments = ([4.,4.,4.,90.,90.,90.], ['x,y,z'], [(1,.25,.25,.25)], [1], 6.)
# Warm up C++ and Python caches before installing one-allocation failure hooks.
cpplib.VoronoiCalculation(*arguments)
failures = 0
for index in range(1, 200):
    _testcapi.set_nomemory(index, index + 1)
    try:
        cpplib.VoronoiCalculation(*arguments)
    except MemoryError:
        failures += 1
    finally:
        _testcapi.remove_mem_hooks()
assert failures > 0
assert cpplib.VoronoiCalculation(*arguments)['voronoi_cells']
print('allocation paths checked')
"""
        result = _run_isolated(script, extension=EXTENSION, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
