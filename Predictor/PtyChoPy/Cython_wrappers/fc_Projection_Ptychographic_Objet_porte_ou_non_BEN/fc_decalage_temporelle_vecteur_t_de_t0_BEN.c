/*********************************************************************************************************
* Date: 6 juillet 2018																					 *
* Auteur: Benoit Brizard																				 *
*																										 *
* Fonction cr�ant une matrice M_decalage � partir du vecteur V defini sur t,  d�cal�e de (t0*sign)		 *
* � chaque ligne. On note que sign ne peut prendre que 2 valeurs->sign = +-1.							 *
* IMPORTANT : Il faut que le pas de t0 soit multiple de celui de t !!									 *
* Il s'agit vraiment du m�me principe que dans fc_decalage_temporelle_Matrice_t_de_t0_BEN, sauf qu'ici   *
* on n'a qu'un vecteur de donn� -> plus simple														     *
* INPUTS : t0, V_t, t, sign																				 *
* OUTPUTS :  M_decalage_t0																			     *
*********************************************************************************************************/
#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
// #include <mex.h>
// #include <Python.h>

// Your function definitions

/*********************************************************************************************************
* Strat�gie : � la place de vraiment ajouter des colonnes de part et d'autres d'une matrice comme dans	 *
* l'algorithme Matlab, on va indexer de mani�re plus efficace. On ne cr�era pas une �norme matrice pour  *
* n'en s�lectionner qu'une partie � la fin : on calcule seulement les indices qu'on garde � la fin, ce qui
* �vite de faire des calculs inutiles. � la fin du script Matlab, on s�lectionne le output ainsi :		 *
* M_decalage_t0 = M_decal_loc(:, N_pts_plus/2+1:1:N+N_pts_plus/2); . Ce que ce programme fait est de	 *
* v�rifier si nos indices se trouve bien dans [N_pts_plus/2+1:1:N+N_pts_plus/2] et si oui de faire le	 *
* calcul appropri�.																						 *
**********************************************************************************************************/
void fc_decalage_temporelle_vecteur_t_de_t0_BEN(double dt0, double dt, int N_pts_plus, double* V_t, double* IM_V_t, double*sign, int N0, int N, double*M_decalage_t0, double*IM_M_decalage_t0)
{
	int ind;
	int moitie = N_pts_plus / 2;
	if (IM_V_t == NULL)
	{
		
		if (*sign == 1.0)
		{
			for (int row = 0; row < N0; row++)
			{
				ind = (int)(row * dt0 / dt) + 1;
				for (int col = 0; col < ind - 1; col++)
				{
					if (col >= moitie && col <= moitie + N - 1)
						M_decalage_t0[row + N0 * (col - moitie)] = V_t[0];

				}
				for (int n = ind - 1; n < ind + N - 2; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
						M_decalage_t0[row + N0 * (n - moitie)] = V_t[n - ind + 1];

				}
				for (int n = ind + N - 2; n < N + N_pts_plus; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
						M_decalage_t0[row + N0 * (n - moitie)] = V_t[N - 1];

				}

			}
		}
		else if (*sign == -1.0)
		{
			for (int row = 0; row < N0; row++)
			{
				ind = N + N_pts_plus - row * (int)(dt0 / dt);
				for (int col = 0; col < ind - N; col++)
				{
					if (col >= moitie && col <= moitie + N - 1)
						M_decalage_t0[row + N0 * (col - moitie)] = V_t[0];

				}
				for (int n = ind - N; n < ind - 1; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
						M_decalage_t0[row + N0 * (n - moitie)] = V_t[n - ind + N];

				}
				for (int n = ind - 1; n < N + N_pts_plus; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
						M_decalage_t0[row + N0 * (n - moitie)] = V_t[N - 1];

				}

			}
		}
		else if (*sign != 1.0 && *sign != -1.0)
		{
		
			// mexErrMsgTxt("Wrong format input 'sign'\n");
		}

	}


	if (IM_V_t != NULL)
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
						M_decalage_t0[row + N0 * (col - moitie)] = V_t[0];
						IM_M_decalage_t0[row + N0 * (col - moitie)] = IM_V_t[0];
					}

				}
				for (int n = ind - 1; n < ind + N - 2; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						M_decalage_t0[row + N0 * (n - moitie)] = V_t[n - ind + 1];
						IM_M_decalage_t0[row + N0 * (n - moitie)] = IM_V_t[n - ind + 1];
					}

				}
				for (int n = ind + N - 2; n < N + N_pts_plus; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						M_decalage_t0[row + N0 * (n - moitie)] = V_t[N - 1];
						IM_M_decalage_t0[row + N0 * (n - moitie)] = IM_V_t[N - 1];
					}

				}

			}
		}

		else if (*sign == -1.0)
		{
			for (int row = 0; row < N0; row++)
			{
				ind = N + N_pts_plus - row * (int)(dt0 / dt);
				for (int col = 0; col < ind - N; col++)
				{
					if (col >= moitie && col <= moitie + N - 1)
					{
						M_decalage_t0[row + N0 * (col - moitie)] = V_t[0];
						IM_M_decalage_t0[row + N0 * (col - moitie)] = IM_V_t[0];
					}

				}
				for (int n = ind - N; n < ind - 1; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						M_decalage_t0[row + N0 * (n - moitie)] = V_t[n - ind + N];
						IM_M_decalage_t0[row + N0 * (n - moitie)] = IM_V_t[n - ind + N];
					}

				}
				for (int n = ind - 1; n < N + N_pts_plus; n++)
				{
					if (n >= moitie && n <= moitie + N - 1)
					{
						M_decalage_t0[row + N0 * (n - moitie)] = V_t[N - 1];
						IM_M_decalage_t0[row + N0 * (n - moitie)] = IM_V_t[N - 1];
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
//Sert � cr�er un .mexw64, ce qui m'aide a debug, mais en temps et lieu je vais appeler fc_decalage_temporelle_vecteur dans fc_contrainte_ptyco directement
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs != 1 || nrhs != 4)
	{
		mexErrMsgTxt("Wrong number of arguments during call");
		return;
	}

	double * IM_V_t;

	if (mxIsComplex(prhs[1]))
	{
	    IM_V_t= mxGetPi(prhs[1]);
	}
	else
	{
		IM_V_t = NULL;
	}


	int m = mxGetM(prhs[1]); //1;
	if (m!=1)
	{
		mexErrMsgTxt("Wrong type of input, entry must me be a row vector");
		return;
	}
	int N0 = mxGetN(prhs[0]);
	int N = mxGetN(prhs[1]); //201;

	
	//INPUTS
	double *V_t = mxGetPr(prhs[1]);
	double*t0 = mxGetPr(prhs[0]);
	double*t = mxGetPr(prhs[2]);
	double dt, dt0;
	dt = t[9] - t[8]; dt0 = t0[9] - t0[8];
	int N_pts_plus = (int)((N0 - 1)* dt0 / dt);

	double* sign = mxGetPr(prhs[3]);

	//OUTPUTS
	double *M_decalage_t0;
	double* IM_M_decalage_t0;
	if (IM_V_t == NULL)
	{
		plhs[0] = mxCreateDoubleMatrix(N0, N, mxREAL);
		M_decalage_t0 = mxGetPr(plhs[0]);
		IM_M_decalage_t0 = NULL;
	}
	else
	{
		plhs[0] = mxCreateDoubleMatrix(N0, N, mxCOMPLEX);
		M_decalage_t0 = mxGetPr(plhs[0]);
		IM_M_decalage_t0 = mxGetPi(plhs[0]);
	}

	fc_decalage_temporelle_vecteur_t_de_t0_BEN(dt0, dt, N_pts_plus, V_t,IM_V_t, sign, N0, N, M_decalage_t0, IM_M_decalage_t0);
	return;
}

*/