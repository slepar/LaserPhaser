#ifndef FC_DECALAGE_TEMPORELLE_VECTEUR_T_DE_T0_BEN_H
#define FC_DECALAGE_TEMPORELLE_VECTEUR_T_DE_T0_BEN_H

void fc_decalage_temporelle_vecteur_t_de_t0_BEN(double dt0, double dt, int N_pts_plus, 
	double* V_t, double* IM_V_t, double*sign, int N0, int N, double*M_decalage_t0,
	double*IM_M_decalage_t0);

#endif