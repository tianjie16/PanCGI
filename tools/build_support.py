from pathlib import Path
import shutil

from setuptools.command.build_py import build_py


class BuildPy(build_py):
    def run(self):
        super().run()
        root = Path(__file__).resolve().parents[1]
        target = Path(self.build_lib) / "pancgi_app"
        for name in ("README.md", "DEPENDENCIES.md", "LICENSE", "NOTICE", "environment.yml"):
            shutil.copy2(root / name, target / name)
        for name in ("docs", "examples", "assets", "resources"):
            shutil.copytree(root / name, target / name, dirs_exist_ok=True)
