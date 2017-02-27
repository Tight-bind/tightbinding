from distutils.core import setup
from Cython.Build import cythonize
#from distutils.extension import Extension

#ext_modules = [
#        Extension("Bakerbind",
#                  sources=["get_input_files.pyx",
#                           "HamiltonianMatrix.pyx"],
#                  libraries=["m"]
#        )
#]

ext_modules = ["get_input_files.pyx",
               "HamiltonianMatrix.pyx",
               "slater_koster.pyx"]
setup(
     ext_modules=cythonize(ext_modules)
)
