import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest


EXTENSION = pathlib.Path(os.environ["CPPLIB_TEST_EXTENSION"]) if "CPPLIB_TEST_EXTENSION" in os.environ else (
    pathlib.Path(__file__).resolve().parents[3]
    / "module"
    / ("cpplib.pyd" if sys.platform == "win32" else "cpplib.so")
)

BONDS = "1\n1 1 0.65 1.2\n"
ATOMS = "atoms = [(1, 0.0, 0.0, 0.0), (1, 0.74, 0.0, 0.0)]"


def _run_isolated(script, *, extension, cwd=None, timeout=30):
    """Run the script in a fresh interpreter.

    Each singleton holds per-process state, so observing "the very first
    call" requires an isolated process.
    """
    command = [
        sys.executable, "-I", "-c",
        textwrap.dedent(script),
        str(extension),
    ]
    return subprocess.run(
        command, cwd=cwd, capture_output=True, text=True, timeout=timeout,
    )


# Shared prelude. Only imports the extension — nothing else. The binding
# resolves data next to the .pyd/.so itself: on the first function call it
# reads `__file__`, takes its parent directory, and hands that path to
# Distances::instance / ElectronDensity::instance.
_LOAD_EXTENSION = """
import importlib.util
import os
import pathlib
import sys

ext = pathlib.Path(sys.argv[1]).resolve()

if sys.platform == "win32":
    dll_directories = [os.add_dll_directory(str(ext.parent))]
    if os.environ.get("CPPLIB_TEST_DLL_DIR"):
        dll_directories.append(os.add_dll_directory(os.environ["CPPLIB_TEST_DLL_DIR"]))

spec = importlib.util.spec_from_file_location("cpplib", ext)
cpplib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cpplib)

print("import OK", flush=True)
"""


class DistanceLoadingTests(unittest.TestCase):
    """Laziness of cpplib singletons.

    1. Original directory: import and calls succeed.
    2. Extension copied WITHOUT its data: import must not touch the
       filesystem and must not raise.
    3. Same copy, function is called: missing data must surface as
       RuntimeError containing "FILE_NOT_FOUND".
    4. Copy together with its data: works.
    """

    @classmethod
    def setUpClass(cls):
        if not EXTENSION.is_file():
            raise unittest.SkipTest(f"Extension not found: {EXTENSION}")

    def _copy_extension(self, destination):
        """Copy the extension into `destination`; return its full path."""
        destination.mkdir(parents=True, exist_ok=True)
        target = destination / ("cpplib" + EXTENSION.suffix)
        shutil.copy2(EXTENSION, target)
        return target

    # ---------- 1. Run from the original directory ---------------------

    def test_original_directory_import_and_call(self):
        script = _LOAD_EXTENSION + f"""
{ATOMS}
assert len(cpplib.GenBonds(atoms)["bonds"]) == 1
assert len(cpplib.GenBondsEx(atoms)["bonds"]) == 1
print("OK")
"""
        result = _run_isolated(
            script, extension=EXTENSION, cwd=EXTENSION.parent,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("import OK", result.stdout)
        self.assertIn("OK", result.stdout)

    # ---------- 2. Remote directory, NO data: import must not raise ----

    def test_import_from_remote_directory_does_not_raise(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            extension = self._copy_extension(root / "remote-\u0434\u0430\u043d\u043d\u044b\u0435")

            cwd = root / "cwd"
            cwd.mkdir()

            result = _run_isolated(_LOAD_EXTENSION, extension=extension, cwd=cwd)
            self.assertEqual(
                result.returncode, 0, result.stdout + result.stderr,
            )
            self.assertIn("import OK", result.stdout)
            self.assertNotIn(
                "FILE_NOT_FOUND", result.stdout + result.stderr,
            )

    # ---------- 3. Remote directory, NO data: call must fail -----------

    def test_call_after_copy_without_data_raises_file_not_found(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            extension = self._copy_extension(root / "remote-\u0434\u0430\u043d\u043d\u044b\u0435")

            cwd = root / "cwd"
            cwd.mkdir()

            script = _LOAD_EXTENSION + f"""
{ATOMS}
for method in (cpplib.GenBonds, cpplib.GenBondsEx):
    try:
        method(atoms)
    except RuntimeError as error:
        assert "FILE_NOT_FOUND" in str(error), f"Unexpected: {{error}}"
    else:
        raise AssertionError(
            f"{{method.__name__}} did not raise without data files"
        )
print("OK")
"""
            result = _run_isolated(script, extension=extension, cwd=cwd)
            self.assertEqual(
                result.returncode, 0, result.stdout + result.stderr,
            )
            self.assertIn("OK", result.stdout)

    # ---------- 4. Copy WITH data works --------------------------------

    def test_relocated_extension_with_data_works(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            extension = self._copy_extension(root / "remote-\u0434\u0430\u043d\u043d\u044b\u0435")
            extension_dir = extension.parent
            (extension_dir / "BondLength.ini").write_text(
                BONDS, encoding="utf-8",
            )

            script = _LOAD_EXTENSION + f"""
{ATOMS}
assert len(cpplib.GenBonds(atoms)["bonds"]) == 1
assert len(cpplib.GenBondsEx(atoms)["bonds"]) == 1
print("OK")
"""
            result = _run_isolated(
                script, extension=extension, cwd=extension_dir,
            )
            self.assertEqual(
                result.returncode, 0, result.stdout + result.stderr,
            )
            self.assertIn("OK", result.stdout)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        EXTENSION = pathlib.Path(sys.argv.pop(1)).resolve()
    unittest.main()