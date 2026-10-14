"""Run with: python test_distance_loading.py /absolute/path/to/cpplib.so (or .pyd)."""

import pathlib
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest


EXTENSION = None
BONDS = "1\n1 1 0.65 1.2\n"


class DistanceLoadingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if EXTENSION is None:
            raise unittest.SkipTest("Pass a built cpplib extension when running this script directly")

    def run_extension(self, script, data=None):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            module_dir = root / "relocated" / "module"
            module_dir.mkdir(parents=True)
            extension = module_dir / ("cpplib" + EXTENSION.suffix)
            shutil.copy2(EXTENSION, extension)
            if data is not None:
                (module_dir / "BondLength.ini").write_text(data)
            cwd = root / "working"
            (cwd / "module").mkdir(parents=True)
            # Valid decoys must never mask a missing file beside the extension.
            (cwd / "BondLength.ini").write_text(BONDS)
            (cwd / "module" / "BondLength.ini").write_text(BONDS)
            setup = """
                import importlib.util
                import pathlib
                import sys
                module_file = pathlib.Path(sys.argv[1])
                spec = importlib.util.spec_from_file_location("cpplib", module_file)
                cpplib = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(cpplib)
                data_file = module_file.with_name("BondLength.ini")
                atoms = [(1, 0.0, 0.0, 0.0), (1, 0.74, 0.0, 0.0)]
            """
            result = subprocess.run(
                [sys.executable, "-I", "-c", textwrap.dedent(setup) + textwrap.dedent(script), str(extension)],
                cwd=cwd, capture_output=True, text=True, timeout=30,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_relocated_extension_loads_bonds(self):
        self.run_extension("""
            assert len(cpplib.GenBonds(atoms)["bonds"]) == 1
            assert len(cpplib.GenBondsEx(atoms)["bonds"]) == 1
        """, BONDS)

    def test_missing_data_raises_and_can_be_retried(self):
        self.run_extension("""
            for method in (cpplib.GenBonds, cpplib.GenBondsEx):
                try:
                    method(atoms)
                except RuntimeError as error:
                    assert "FILE_NOT_FOUND" in str(error), str(error)
                else:
                    raise AssertionError("Missing data did not raise")
            data_file.write_text("1\\n1 1 0.65 1.2\\n")
            assert len(cpplib.GenBonds(atoms)["bonds"]) == 1
        """)

    def test_malformed_data_raises_until_repaired(self):
        self.run_extension("""
            for data in ("-1\\n", "1\\n1 1 0.65", "1\\n0 1 0.65 1.2\\n"):
                data_file.write_text(data)
                try:
                    cpplib.GenBonds(atoms)
                except RuntimeError:
                    pass
                else:
                    raise AssertionError("Malformed data did not raise")
            data_file.write_text("1\\n1 1 0.65 1.2\\n")
            assert len(cpplib.GenBonds(atoms)["bonds"]) == 1
        """)

    def test_other_cpp_exceptions_become_python_errors(self):
        self.run_extension("""
            try:
                cpplib.Cluster([10., 10., 10., 90., 90., 90.], ["x,y,z"], atoms,
                               [(0., 0., 0., -1.)], 3.)
            except RuntimeError as error:
                assert "Anchor radius must be positive" in str(error), str(error)
            else:
                raise AssertionError("Invalid anchor did not raise")
        """, BONDS)


if __name__ == "__main__":
    EXTENSION = pathlib.Path(sys.argv.pop(1)).resolve()
    unittest.main()
