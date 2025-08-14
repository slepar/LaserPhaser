
from Software.Predictor.PtyChoPy.Cython_wrappers.fc_Projection_Ptychographic_Objet_porte_ou_non_BEN import projection_ptychographic
from Software.Predictor.PtyChoPy.Cython_wrappers.fc_projection_experimentale_SYD import fc_projection_experimentale_SYD

def fc_iteration_HIO(Trace_Int_XP, M_Produit_t_t0_it, t, t0, k_max, Case_objet):

    # 1A) Pi_OS(Prod_i)
    M_Produit_t_t0_it_inter1A, Obj_t_it, Chp_t_it = projection_ptychographic.call_fc_Projection_Ptychographic_Objet_porte_ou_non_BEN(
                t, t0, k_max, Case_objet, M_Produit_t_t0_it.T)   
    M_Produit_t_t0_it_inter1A = M_Produit_t_t0_it_inter1A.T

    # 1B) Pi_XP(Prod_i)
    M_Produit_t_t0_it_inter1B, Trace_amp_omg_t0_it, Err_it = fc_projection_experimentale_SYD.fc_projection_experimentale_py(
        M_Produit_t_t0_it, Trace_Int_XP)

    # 2A) Pi_XP A
    M_Produit_t_t0_it_inter2A, pbll1, pbll2 = fc_projection_experimentale_SYD.fc_projection_experimentale_py(
        M_Produit_t_t0_it_inter1A, Trace_Int_XP)

    # 2B) Pi_OS B
    M_Produit_t_t0_it_inter2B, pbll1, pbll2 = projection_ptychographic.call_fc_Projection_Ptychographic_Objet_porte_ou_non_BEN(
                t, t0, k_max, Case_objet, M_Produit_t_t0_it_inter1B.T)
    M_Produit_t_t0_it_inter2B = M_Produit_t_t0_it_inter2B.T

    # Sum
    M_Produit_t_t0_it_plus2 = M_Produit_t_t0_it + (M_Produit_t_t0_it_inter2A - M_Produit_t_t0_it_inter2B)

    return Trace_amp_omg_t0_it, M_Produit_t_t0_it_plus2, Obj_t_it, Chp_t_it, Err_it