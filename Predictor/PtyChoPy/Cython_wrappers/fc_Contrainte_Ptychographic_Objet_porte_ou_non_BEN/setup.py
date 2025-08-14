from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

# Define extension modules
extensions = [
    Extension(
        "contrainte_ptychographic",
        sources=["contrainte_ptychographic.pyx", "fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN.c",
                 "fc_decalage_temporelle_Matrice_t_de_t0_BEN.c",
                 "fc_decalage_temporelle_vecteur_t_de_t0_BEN.c",
                 "fc_somme_sur_les_dt0_decales_Matrice_de_t_decalee_de_t0_BEN.c"],

        include_dirs=[np.get_include()],  # Include NumPy headers
    )
]

# Setup configuration
setup(
    name="your_package_name",
    ext_modules=cythonize(extensions),
    include_dirs=[np.get_include()],  # Include NumPy headers for all extension modules
)
