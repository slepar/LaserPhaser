from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension(
        "projection_ptychographic", 
        sources=["projection_ptychographic.pyx", 
                 "fc_Projection_Ptychographic_Objet_porte_ou_non_BEN.c",
                 "fc_decalage_temporelle_vecteur_t_de_t0_BEN.c",
                 "fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN.c",
                 "fc_decalage_temporelle_Matrice_t_de_t0_BEN.c",
                 "fc_somme_sur_les_dt0_decales_Matrice_de_t_decalee_de_t0_BEN.c"],  
        include_dirs=[np.get_include()],  # Include NumPy headers
    )
]

setup(
    name="projection_ptychographic",  # Replace "your_package_name" with the desired package name
    ext_modules=cythonize(ext_modules),
    include_dirs=[np.get_include()],  # Include NumPy headers
)
