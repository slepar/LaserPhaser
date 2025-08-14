# Import necessary modules
import numpy as np
cimport numpy as np

# Declare the C function
cdef extern from "fc_decalage_temporelle_vecteur_t_de_t0_BEN.h":
    void fc_decalage_temporelle_vecteur_t_de_t0_BEN(
        double dt0, double dt, int N_pts_plus,
        double* V_t, double* IM_V_t, double* sign,
        int N0, int N,
        double* M_decalage_t0, double* IM_M_decalage_t0)

# Define the Python function that wraps the C function
def call_fc_decalage_temporelle_vecteur_t_de_t0_BEN(
                                np.ndarray[np.float64_t, ndim=1] t, 
                                np.ndarray[np.float64_t, ndim=1] t0, 
                                np.ndarray[np.complex128_t, ndim=2] V_t, 
                                double sign):

    # Calculate dt0, dt, N_pts_plus, N0, and N
    cdef int N0 = int(t0.shape[0])       # t is 1D here so cols is rows
    cdef int N = int(V_t.shape[1])
    cdef double dt0 = t0[9] - t0[8]
    cdef double dt = t[9] - t[8]
    cdef int N_pts_plus = int((N0 - 1) * dt0 / dt)

    # declare output/input vectors
    cdef np.ndarray[np.float64_t, ndim=2] M_decalage_t0
    cdef np.ndarray[np.float64_t, ndim=2] IM_M_decalage_t0
    cdef np.ndarray[np.float64_t, ndim=2] REAL_V_t = V_t.real.copy()
    cdef np.ndarray[np.float64_t, ndim=2] IM_V_t = V_t.imag.copy()

    if np.iscomplexobj(V_t):
        M_decalage_t0 = np.empty((N, N0), dtype=np.float64)
        IM_M_decalage_t0 = np.empty((N, N0), dtype=np.float64)

        fc_decalage_temporelle_vecteur_t_de_t0_BEN(
                                dt0, dt, N_pts_plus,
                                &REAL_V_t[0, 0], &IM_V_t[0, 0], &sign, N0, N,
                                &M_decalage_t0[0, 0], &IM_M_decalage_t0[0, 0])
        return M_decalage_t0 + 1j * IM_M_decalage_t0
        
    else:
        M_decalage_t0 = np.empty((N, N0), dtype=np.float64)
        fc_decalage_temporelle_vecteur_t_de_t0_BEN(
                                            dt0, dt, N_pts_plus,
                                            &REAL_V_t[0, 0], NULL, &sign, N0, N,
                                            &M_decalage_t0[0, 0], NULL)
        return M_decalage_t0