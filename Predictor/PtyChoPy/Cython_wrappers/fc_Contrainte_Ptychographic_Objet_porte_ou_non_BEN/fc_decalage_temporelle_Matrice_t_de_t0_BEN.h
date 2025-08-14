#ifndef FC_DECALAGE_TEMPORELLE_MATRICE_T_DE_T0_BEN_H
#define FC_DECALAGE_TEMPORELLE_MATRICE_T_DE_T0_BEN_H

void fc_decalage_temporelle_Matrice_t_de_t0_BEN(double dt0, double dt, int N_pts_plus, 
	double* REAL_M0, double* IM_M0, double*sign, int N0, int N,
	double*REAL_M_decal_loc, double*IM_M_decal_loc);

#endif
