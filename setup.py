from distutils.core import setup
from Cython.Build import cythonize
import numpy as np
from distutils.extension import Extension

ext_modules = [
        Extension("Get_input",
                  sources=["get_input_files.pyx"],
                  libraries=["m"]
        ),
            Extension("HamiltonianMatrix",
                      sources=["HamiltonianMatrix.pyx"],
                      libraries=["m"],
                      include_dirs=[np.get_include()]
        ),
            Extension("slater_koster",
                          sources=["slater_koster.pyx"],
                          libraries=["m"]
        )

]

#ext_modules = ["get_input_files.pyx",
#               "HamiltonianMatrix.pyx",
#               "slater_koster.pyx"]
setup(
     ext_modules=cythonize(ext_modules)
     , include_dirs=[np.get_include()]
)
