# Import necessary modules
import numpy as np
cimport numpy as np
from libc.math cimport floor

# Declare the C functions
cdef extern from "fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN.h":
    void fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(
        double max_t0, double dt0, double dt, double* t,
        int N_pts_plus_somme, int N_pts_plus_decalage, int N_dt0_s_dt,
        double k_max, double case_object,
        double* REAL_Produit_t_t0, double* IM_Produit_t_t0,
        int N0, int N, double* REAL_Obj_t, double* IM_Obj_t,
        double* REAL_Chp_t, double* IM_Chp_t)

def call_fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(
        np.ndarray[np.float64_t, ndim=1] t, 
        np.ndarray[np.float64_t, ndim=1] t0,
        double k_max, 
        double case_object,
        np.ndarray[np.complex128_t, ndim=2] Produit_t_t0):

    # checking input shapes
    if Produit_t_t0.shape[0] != t.shape[0] or Produit_t_t0.shape[1] != t0.shape[0]:
        raise ValueError(f"Matrix must be ({t.shape[0]}, {t0.shape[0]})")
    
    # calculating unassigned variables
    cdef int N0 = int(Produit_t_t0.shape[1])
    cdef int N = int(Produit_t_t0.shape[0])
    cdef double max_t0 = t0[0]
    cdef int i
    for i in range(1, N0):
        if t0[i] > max_t0:
            max_t0 = t0[i]
    cdef double dt0 = t0[9] - t0[8]
    cdef double dt = t[9] - t[8]
    cdef int N_dt0_s_dt = int(dt0 / dt)          
    cdef int N_pts_plus_somme = int(floor((N_dt0_s_dt + 1.0) / 2.0))
    cdef int N_pts_plus_decalage = int((N0 - 1)* N_dt0_s_dt)
    cdef np.ndarray[np.float64_t, ndim=2] t2D = t.reshape(1, -1)  

    # assigning real and imaginary parts of input
    cdef np.ndarray[np.float64_t, ndim=2] REAL_Produit_t_t0 = Produit_t_t0.real.copy()
    cdef np.ndarray[np.float64_t, ndim=2] IM_Produit_t_t0 = Produit_t_t0.imag.copy() if np.iscomplexobj(Produit_t_t0) else None
    
    # creating output vectors
    cdef np.ndarray[np.complex128_t, ndim=2] Obj_t 
    cdef np.ndarray[np.complex128_t, ndim=2] Chp_t 
    cdef np.ndarray[np.float64_t, ndim=2] REAL_Obj_t 
    cdef np.ndarray[np.float64_t, ndim=2] IM_Obj_t
    cdef np.ndarray[np.float64_t, ndim=2] REAL_Chp_t 
    cdef np.ndarray[np.float64_t, ndim=2] IM_Chp_t
    
    if np.iscomplexobj(Produit_t_t0):
        REAL_Obj_t = np.empty((1, N), dtype=np.float64)
        IM_Obj_t = np.empty((1, N), dtype=np.float64)

        REAL_Chp_t = np.empty((1, N), dtype=np.float64)
        IM_Chp_t = np.empty((1, N), dtype=np.float64)

        fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(
            max_t0, dt0, dt, &t2D[0, 0], 
            N_pts_plus_somme, N_pts_plus_decalage, N_dt0_s_dt, 
            k_max, case_object,
            &REAL_Produit_t_t0[0, 0], &IM_Produit_t_t0[0, 0],
            N0, N, 
            &REAL_Obj_t[0, 0], &IM_Obj_t[0, 0], 
            &REAL_Chp_t[0, 0], &IM_Chp_t[0, 0])

        return REAL_Obj_t + 1j * IM_Obj_t, REAL_Chp_t + 1j * IM_Chp_t

    else:
        REAL_Obj_t = np.empty((1, N), dtype=np.float64)
        REAL_Chp_t = np.empty((1, N), dtype=np.float64)
        fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(
            max_t0, dt0, dt, &t2D[0, 0], 
            N_pts_plus_somme, N_pts_plus_decalage, N_dt0_s_dt, 
            k_max, case_object,
            &REAL_Produit_t_t0[0, 0], &IM_Produit_t_t0[0, 0],
            N0, N, 
            &REAL_Obj_t[0, 0], NULL, 
            &REAL_Chp_t[0, 0], NULL)

        return REAL_Obj_t, REAL_Chp_t