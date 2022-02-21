import sys
from pybind11 import get_cmake_dir
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
__version__ = "0.0.1"

ext_modules = [
    Pybind11Extension(
        "d5lidar",
        ["src/d5lidar.cpp", 
         "src/BinFile.cpp",
         "src/thirdparty/miniz.c"],
        ),
    ]

setup(
    name="d5lidar",
    version=__version__,
    author="M. Grady Saunders",
    author_email="mgscis@rit.edu",
    url="https://github.com/mgradysaunders/d5lidar"
    description="DIRSIG5 Lidar utilities.",
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.6"
)
