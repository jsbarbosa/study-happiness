from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "core",
        ["core.pyx"],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp', '-lgomp'],
    )
]

setup(
    ext_modules = cythonize(ext_modules),
    include_dirs = [numpy.get_include()]
)
