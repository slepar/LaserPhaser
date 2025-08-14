from Software.Predictor.PtyChoPy.Cython_wrappers.fc_Projection_Ptychographic_Objet_porte_ou_non_BEN import projection_ptychographic
from Software.Predictor.PtyChoPy.Cython_wrappers.fc_projection_experimentale_SYD import fc_projection_experimentale_SYD

def fc_iteration_ER(Trace_Int_XP, M_Produit_t_t0_it, t, t0, k_max, Case_objet):

    # 1) Pi_XP : Prod_i -> Prod_i+1
    M_Produit_t_t0_it_plus_1, Trace_amp_omg_t0_it, Err_it = fc_projection_experimentale_SYD.fc_projection_experimentale_py(
        M_Produit_t_t0_it, Trace_Int_XP)

    # 2) Pi_OS : Prod_i+1 -> Prod_i+2
    M_Produit_t_t0_it_plus2, Obj_t_it_p2, Chp_t_it_2 = projection_ptychographic.call_fc_Projection_Ptychographic_Objet_porte_ou_non_BEN(
        t, t0, k_max, Case_objet, M_Produit_t_t0_it_plus_1.T)
    M_Produit_t_t0_it_plus2 = M_Produit_t_t0_it_plus2.T

    return Trace_amp_omg_t0_it, M_Produit_t_t0_it_plus2, Obj_t_it_p2, Chp_t_it_2, Err_it

