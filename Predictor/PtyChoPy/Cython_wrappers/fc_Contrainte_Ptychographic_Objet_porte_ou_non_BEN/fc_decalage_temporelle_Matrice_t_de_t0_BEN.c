/*********************************************************************************************************
* Date: 6 juillet 2018																					 *
* Auteur: Benoit Brizard																				 *
*																										 *
* Fonction cr�ant une matrice M_decalage � partir d'une matrice M definie sur t et t0 d�cal�e de(t0*sign)*  
* � chaque ligne. On note que sign ne peut prendre que 2 valeurs->sign = +-1.							 *
* IMPORTANT : Il faut que le pas de t0 soit multiple de celui de t !!									 *
* INPUTS :  t0, M0, t, sign																			     *
* OUTPUTS : Produit_t_t0_decal																		     *
*********************************************************************************************************/


#define _CRT_SECURE_NO_WARNINGS
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <stdlib.h>
#include <stdio.h>
// #include <mex.h>
#include "fc_decalage_temporelle_Matrice_t_de_t0_BEN.h"


/*********************************************************************************************************
* Strat�gie : � la place de vraiment ajouter des colonnes de part et d'autres d'une matrice comme dans	 *
* l'algorithme Matlab, on va indexer de mani�re plus efficace. On ne cr�era pas une �norme matrice pour  *
* n'en s�lectionner qu'une partie � la fin : on calcule seulement les indices qu'on garde � la fin, ce qui
* �vite de faire des calculs inutiles. � la fin du script Matlab, on s�lectionne le output ainsi :		 *
* M_decalage_t0 = M_decal_loc(:, N_pts_plus/2+1:1:N+N_pts_plus/2); . Ce que ce programme fait est de	 *
* v�rifier si nos indices se trouve bien dans [N_pts_plus/2+1:1:N+N_pts_plus/2] et si oui de faire le	 *
* calcul appropri�.																						 *
**********************************************************************************************************/


void fc_decalage_temporelle_Matrice_t_de_t0_BEN(double dt0,double dt,int N_pts_plus, double* REAL_M0, double* IM_M0, 
												double*sign, int N0, int N,
												double*REAL_Produit_t_t0_decal, double*IM_Produit_t_t0_decal)
{
	int ind; int moitie = N_pts_plus / 2;
	if (IM_M0 != NULL)
	{
		if (*sign == 1.0)
		{
			for (int row = 0; row < N0; row++)
			{
				ind = (int)(row * dt0 / dt) + 1;
				for (int col = 0; col < ind - 1; col++)
				{
					if (col >= moitie && col <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (col-moitie)] = REAL_M0[row];
						IM_Produit_t_t0_decal[row + N0 * (col - moitie)] = IM_M0[row];
					}
				}
				for (int n = ind - 1; n < ind + N - 2; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (n-moitie)] = REAL_M0[row + N0 *(n - ind + 1)];
						IM_Produit_t_t0_decal[row + N0 * (n - moitie)] = IM_M0[row + N0 * (n - ind + 1)];
					}
				}
				for (int n = ind + N - 2; n < N + N_pts_plus; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (n - moitie)] = REAL_M0[row + N0 * (N - 1)];
						IM_Produit_t_t0_decal[row + N0 * (n - moitie)] = IM_M0[row + N0 * (N - 1)];
					}
				}

			}

		}

		else if (*sign == -1.0)
		{
			for (int row = 0; row < N0; row++)
			{
				ind = N + N_pts_plus - row * ((int)(dt0 / dt));
				for (int col = 0; col < ind - N; col++)
				{
					if (col >= moitie && col <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (col - moitie)] = REAL_M0[row];
						IM_Produit_t_t0_decal[row + N0 * (col - moitie)] = IM_M0[row];
					}
				}
				for (int n = ind - N; n < ind - 1; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (n - moitie)] = REAL_M0[row + N0 * (n - ind + N)];
						IM_Produit_t_t0_decal[row + N0 * (n - moitie)] = IM_M0[row + N0 * (n - ind + N)];
					}
				}
				for (int n = ind - 1; n < N + N_pts_plus; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (n - moitie)] = REAL_M0[row + N0 * (N - 1)];
						IM_Produit_t_t0_decal[row + N0 * (n - moitie)] = IM_M0[row + N0 * (N - 1)];
					}
				}

			}
		}

		else if (*sign != 1.0 && *sign != -1.0)
		{
			// mexErrMsgTxt("Wrong format input 'sign'\n");
		}
	}

	if (IM_M0 == NULL)
	{
		if (*sign == 1.0)
		{
			for (int row = 0; row < N0; row++)
			{
				ind = (int)(row * dt0 / dt) + 1;
				for (int col = 0; col < ind - 1; col++)
				{
					if (col >= moitie && col <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (col - moitie)] = REAL_M0[row];
					}
				}
				for (int n = ind - 1; n < ind + N - 2; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (n - moitie)] = REAL_M0[row + N0 *(n - ind + 1)];
					}
				}
				for (int n = ind + N - 2; n < N + N_pts_plus; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (n - moitie)] = REAL_M0[row + N0 * (N - 1)];
					}
				}

			}

		}

		else if (*sign == -1.0)
		{
			for (int row = 0; row < N0; row++)
			{
				ind = N + N_pts_plus - row * ((int)(dt0 / dt));
				for (int col = 0; col < ind - N; col++)
				{
					if (col >= moitie && col <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (col - moitie)] = REAL_M0[row];
					}
				}
				for (int n = ind - N; n < ind - 1; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (n - moitie)] = REAL_M0[row + N0 * (n - ind + N)];
					}
				}
				for (int n = ind - 1; n < N + N_pts_plus; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						REAL_Produit_t_t0_decal[row + N0 * (n - moitie)] = REAL_M0[row + N0 * (N - 1)];
					}
				}

			}
		}

		else if (*sign != 1.0 && *sign != -1.0)
		{
			// mexErrMsgTxt("Wrong format input 'sign'\n");
		}
	}


	return;
}
/*
//Sert � cr�er un .mexw64, ce qui m'aide a debug, mais en temps et lieu je vais appeler fc_decalage_temporelle_Matrice dans fc_contrainte_ptyco directement
	void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{
		if (nlhs != 1 || nrhs != 4)
		{
			mexErrMsgTxt("Wrong number of arguments during call");
			return;
		}

		double *IM_Produit_t_t0;
		if (!(mxIsComplex(prhs[1])))
		{
			IM_Produit_t_t0 = NULL;
		}
		else
		{
			IM_Produit_t_t0 = mxGetPi(prhs[1]);
		}


		int N0 = mxGetM(prhs[1]); //207;
		int N = mxGetN(prhs[1]); //201;

	
									*
									double *Ar = mxGetPr(A IN); // Real data
									double *Ai = mxGetPi(A IN); // Imaginary data
									And Ar[m + M * n] and Ai[m + M * n] are the real and imaginary parts of A(m + 1, n + 1).
									To create a 2 - D complex array, use
									B OUT = mxCreateDoubleMatrix(M, N, mxCOMPLEX);*/

	/*
		//INPUTS
		double *REAL_Produit_t_t0 = mxGetPr(prhs[1]);	
		double*t0= mxGetPr(prhs[0]);
		double*t = mxGetPr(prhs[2]);
		double dt, dt0;
		dt = t[9] - t[8]; dt0 = t0[9] - t0[8];
		int N_pts_plus = (int) ((N0- 1)* dt0 / dt);
		double* sign = mxGetPr(prhs[3]);

		//OUTPUTS

		double *REAL_Produit_t_t0_decal;
		double *IM_Produit_t_t0_decal;
		if (IM_Produit_t_t0 != NULL)
		{
			plhs[0] = mxCreateDoubleMatrix(N0, N, mxCOMPLEX);
			REAL_Produit_t_t0_decal = mxGetPr(plhs[0]);
			IM_Produit_t_t0_decal = mxGetPi(plhs[0]);
		}
		else
		{
			plhs[0] = mxCreateDoubleMatrix(N0, N, mxREAL);
			REAL_Produit_t_t0_decal = mxGetPr(plhs[0]);
			IM_Produit_t_t0_decal = NULL;
		}
		fc_decalage_temporelle_Matrice_t_de_t0_BEN(dt0, dt, N_pts_plus, REAL_Produit_t_t0, IM_Produit_t_t0, sign, N0, N,
			REAL_Produit_t_t0_decal, IM_Produit_t_t0_decal);
		return;
	}
	*/