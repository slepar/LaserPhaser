import numpy as np
from decalage_temporelle_vecteur import call_fc_decalage_temporelle_vecteur_t_de_t0_BEN

# Example input values
t0 = np.array([0.0, 0.1, 0.2, 0.3], dtype=np.double)
V_t = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.double)
t = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1], dtype=np.double)
sign = 1.0

# Calculating other variables
dt0 = t0[1] - t0[0]  # Calculate dt0
dt = t[1] - t[0]  # Calculate dt
N_pts_plus = (len(t0) - 1) * int(dt0 / dt)
N0 = len(t0)
N = len(V_t)

# Allocate output arrays
M_decalage_t0 = np.zeros((N0, N), dtype=np.double)
IM_M_decalage_t0 = np.zeros((N0, N), dtype=np.double)

# Call the Python function that wraps your C function
M_decalage_t0, IM_M_decalage_t0 = call_fc_decalage_temporelle_vecteur_t_de_t0_BEN(
    dt0, dt, N_pts_plus, V_t, None, sign, N0, N, M_decalage_t0, IM_M_decalage_t0)

print("M_decalage_t0:", M_decalage_t0)
print("IM_M_decalage_t0:", IM_M_decalage_t0)
