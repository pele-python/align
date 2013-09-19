import os
import numpy as np
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension

# # Numpy header files
numpy_lib = os.path.split(np.__file__)[0]
numpy_include = os.path.join(numpy_lib, 'core/include')

setup(
    name="align",
    version='0.0.1',
    description="tools for aligning two molecular structures",
    packages=["align",
              "align._src",
              "align.tests",
             ],
    ext_modules=[
            Extension("align._src.minperm", ["align/_src/minperm.f90"],
                      extra_compile_args=['-Wextra', '-pedantic', '-funroll-loops', '-O2', ],
                     ),
            Extension("align._src._cost_matrix", ["align/_src/_cost_matrix.c"],
                      include_dirs=[numpy_include],
                      extra_compile_args = ['-Wextra','-pedantic','-funroll-loops','-O2',],
                     ),

                 ]
      )
