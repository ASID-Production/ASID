import pathlib
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest


EXTENSION = (
    pathlib.Path(__file__).resolve().parents[3]
    / "module"
    / ("cpplib.pyd" if sys.platform == "win32" else "cpplib.so")
)

BONDS = "1\n1 1 0.65 1.2\n"


class DistanceLoadingTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not EXTENSION.is_file():
            raise unittest.SkipTest(
                f"Extension not found: {EXTENSION}"
            )

    def run_extension(self, script, data=None):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)

            # Каталог, в который переносим расширение
            module_dir = root / "relocated" / "module"
            module_dir.mkdir(parents=True)

            extension = module_dir / ("cpplib" + EXTENSION.suffix)
            shutil.copy2(EXTENSION, extension)

            # Файл данных рядом с расширением, если он нужен тесту
            if data is not None:
                (module_dir / "BondLength.ini").write_text(
                    data, encoding="utf-8"
                )

            # Рабочий каталог с файлами-приманками
            cwd = root / "working"
            (cwd / "module").mkdir(parents=True)

            (cwd / "BondLength.ini").write_text(
                BONDS, encoding="utf-8"
            )
            (cwd / "module" / "BondLength.ini").write_text(
                BONDS, encoding="utf-8"
            )
            (cwd / "ED.txt").write_text(
                "invalid decoy\n", encoding="utf-8"
            )

            # Этот код выполняется в отдельном Python-процессе
            setup = setup = """
                import os
                import pathlib
                import sys
                import importlib.util
                import traceback

                module_file = pathlib.Path(sys.argv[1]).resolve()

                print("Python:", sys.executable)
                print("Extension:", module_file)
                print("Exists:", module_file.exists())
                print("Size:", module_file.stat().st_size)

                dll_directory = os.add_dll_directory(str(module_file.parent))

                try:
                    spec = importlib.util.spec_from_file_location(
                        "cpplib", module_file
                    )
                    cpplib = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(cpplib)
                    print("Extension loaded successfully")
                except Exception:
                    traceback.print_exc()
                    raise
            """

            result = subprocess.run(
                [
                    sys.executable,
                    "-I",
                    "-c",
                    textwrap.dedent(setup) + textwrap.dedent(script),
                    str(extension),
                ],
                cwd=cwd,
                capture_output=True,
                text=True,
                timeout=30,
            )

            self.assertEqual(
                result.returncode,
                0,
                result.stdout + result.stderr,
            )

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

            data_file.write_text(
                "1\\n1 1 0.65 1.2\\n", encoding="utf-8"
            )

            assert len(cpplib.GenBonds(atoms)["bonds"]) == 1
        """)

    def test_ed_is_loaded_lazily_from_the_extension_directory(self):
        self.run_extension("""
            try:
                cpplib.FindCP(
                    [10., 10., 10., 90., 90., 90.],
                    ["x,y,z"],
                    [(1, 0., 0., 0.)],
                    3.
                )
            except RuntimeError as error:
                assert "FILE_NOT_FOUND" in str(error), str(error)
            else:
                raise AssertionError("Missing ED data did not raise")
        """)

    def test_malformed_data_raises_until_repaired(self):
        self.run_extension("""
            for data in (
                "-1\\n",
                "1\\n1 1 0.65",
                "1\\n0 1 0.65 1.2\\n"
            ):
                data_file.write_text(data, encoding="utf-8")

                try:
                    cpplib.GenBonds(atoms)
                except RuntimeError:
                    pass
                else:
                    raise AssertionError("Malformed data did not raise")

            data_file.write_text(
                "1\\n1 1 0.65 1.2\\n", encoding="utf-8"
            )

            assert len(cpplib.GenBonds(atoms)["bonds"]) == 1
        """)

    def test_other_cpp_exceptions_become_python_errors(self):
        self.run_extension("""
            try:
                cpplib.Cluster(
                    [10., 10., 10., 90., 90., 90.],
                    ["x,y,z"],
                    atoms,
                    [(0., 0., 0., -1.)],
                    3.
                )
            except RuntimeError as error:
                assert "Anchor radius must be positive" in str(error)
            else:
                raise AssertionError("Invalid anchor did not raise")
        """ , BONDS)
    def test_extension_loads_from_original_directory(self):
        result = subprocess.run(
            [
                sys.executable,
                "-I",
                "-c",
                """
import importlib.util
import os
import pathlib
import sys

p = pathlib.Path(sys.argv[1]).resolve()

dll_directory = None
if sys.platform == "win32":
    dll_directory = os.add_dll_directory(str(p.parent))

spec = importlib.util.spec_from_file_location("cpplib", p)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

print("OK:", p)
""",
                str(EXTENSION),
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )

        self.assertEqual(
            result.returncode, 0,
            result.stdout + result.stderr
        )

if __name__ == "__main__":
    if len(sys.argv) > 1:
        EXTENSION = pathlib.Path(sys.argv.pop(1)).resolve()

    unittest.main()