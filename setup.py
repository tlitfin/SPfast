from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.0.1"

ext_modules = [
    Pybind11Extension(
        "SPlib",
        ["src/protein2.cc", "src/align1.cc", "src/quatfit_theo.cc",\
            "src/SPfast.cc", "src/quatfit.cc", "src/SPlib.cpp"],
        extra_compile_args=["-std=gnu++11", "-O3", "-shared", "-fopenmp", "-DMP", "-fPIC"],
        extra_link_args=["-std=gnu++11", "-O3", "-shared", "-fopenmp", "-DMP", "-fPIC"]
    ),
]

setup(
    name="SPlib",
    version=__version__,
    author="Thomas Litfin",
    author_email="tomlitfin@gmail.com",
    description="Python bindings for SPfast",
    long_description="",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
