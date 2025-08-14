from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        name="decalage_temporelle_vecteur",
        sources=["decalage_temporelle_vecteur.pyx", "fc_decalage_temporelle_vecteur_t_de_t0_BEN.c"],
        include_dirs=[np.get_include()],  # Include NumPy headers
    )
]

setup(
    ext_modules=cythonize(extensions),
    include_dirs=[np.get_include()]
)
