
#ifndef FC_CONTRAINTE_PTYCHOGRAPHIC_OBJET_PORTE_OU_NON_BEN_H
#define FC_CONTRAINTE_PTYCHOGRAPHIC_OBJET_PORTE_OU_NON_BEN_H

void fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(double max_t0,double dt0, double dt,double * t, int N_pts_plus_somme, 
	int N_pts_plus_decalage, int N_dt0_s_dt,
	double k_max, double case_object,	double *REAL_Produit_t_t0, double *IM_Produit_t_t0, int N0, int N, double *REAL_Obj_t, 
	double *IM_Obj_t, double *REAL_Chp_t, double * IM_Chp_t);

void sum_colonne_normale(double*REAL_Obj_t, double*IM_Obj_t, double*M_Obj_decalage, double*IM_M_Obj_decalage, 
	double*REAL_Produit_t_t0_decal,
	double*IM_Produit_t_t0_decal, int N0, int N);

void mean_abs_filtre(double*r_input, double*im_input, double*V_filtr, int N, double version);

#endif // MY_C_FUNCTIONS_H
