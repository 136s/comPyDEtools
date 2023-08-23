#!/usr/bin/env python

from pathlib import Path
from setuptools import setup, find_packages, Command
import shutil


NAME = "compydetools"
HERE = Path(__file__).resolve().parent


class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""

    CLEAN_FILES = [
        "build",
        "dist",
        "*.pyc",
        "__pycache__",
        "*.egg-info",
        ".eggs",
        "*.egg",
    ]

    user_options = []

    def initialize_options(self) -> None:
        pass

    def finalize_options(self) -> None:
        pass

    def run(self) -> None:
        for pattern in self.CLEAN_FILES:
            for file_path in HERE.glob(pattern):
                if file_path.is_file() or file_path.is_symlink():
                    file_path.unlink()
                elif file_path.is_dir():
                    shutil.rmtree(file_path)


about = {}
with open(HERE.joinpath(NAME, "__version__.py")) as f:
    exec(f.read(), about)

setup(
    name=NAME,
    version=about.get("__version__"),
    description="A Python implementation of a part of compareDEtools.",
    packages=find_packages(),
    package_data={NAME: ["data/*.csv", "data/*.yaml"]},
    cmdclass={"clean": CleanCommand},
)
