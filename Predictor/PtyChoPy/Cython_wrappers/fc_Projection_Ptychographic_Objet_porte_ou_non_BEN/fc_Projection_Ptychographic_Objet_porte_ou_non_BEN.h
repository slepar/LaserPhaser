#ifndef FC_PROJECTION_PTYCHOGRAPHIC_OBJET_PORTE_OU_NON_BEN_H
#define FC_PROJECTION_PTYCHOGRAPHIC_OBJET_PORTE_OU_NON_BEN_H

void circshift(double *REAL_out, double *IM_out, double *REAL_in, double *IM_in, int xdim, int ydim, int xshift);


void fc_Projection_Ptychographic_Objet_porte_ou_non_BEN(double max_t0, double dt0, double dt, double * t, int N_pts_plus_somme, int N_pts_plus_decalage, int N_dt0_s_dt,
	double k_max, double case_object, double *REAL_Produit_t_t0, double *IM_Produit_t_t0, int N0, int N, double *REAL_Obj_t, double *IM_Obj_t, double *REAL_Chp_t,
	double * IM_Chp_t, double*REAL_M_Produit_t_t0_it_plus2, double*IM_M_Produit_t_t0_it_plus2);

#endif 
