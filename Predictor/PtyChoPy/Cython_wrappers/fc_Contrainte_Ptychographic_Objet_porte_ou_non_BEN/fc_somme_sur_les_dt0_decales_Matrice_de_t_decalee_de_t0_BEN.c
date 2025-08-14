/*********************************************************************************************************
* Date: 6 juillet 2018																					 *
* Auteur: Benoit Brizard																				 *
*																										 *
* Fonction cr�ant une matrice M_decalage � partir d'une somme de diff�rentes matrices M_t_t0, qui sont   *
* sur le pas de temps t, d�cal�es de t0.																 *
* INPUTS  : t0, M_t_t0, t 																				 *
* OUTPUTS : M_decalage_t0																				 *
*********************************************************************************************************/


#define _CRT_SECURE_NO_WARNINGS
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
// #include <mex.h>

void fc_somme_sur_les_dt0_decales_Matrice_de_t_decalee_de_t0_BEN(int N_pts_plus, int N_dt0_s_dt, double*M_t_t0, double*IM_M_t_t0,
	double*M_decalage_t0,double*IM_M_decalage_t0, int N0, int N)
{
	
	//Si aucune colonne � ajouter, aucun travail n'est � faire
	if (N_pts_plus == 0)
	{
		return;
	}

	//variables dec : servent pour les index dans les matrices
	int dec,dec1;

	if (IM_M_t_t0 == NULL)
	{
		double * M_t_t0_loc;
		if (N_dt0_s_dt % 2)
		{
			//ajout de N-1 colonnes � gauche et N colonnes � droite de la matrice donn�e en input
			M_t_t0_loc = malloc(sizeof(double)*N0*(N + 2 * N_pts_plus - 1));

			//cr�er d'abord M_t_t0_loc en ajoutant les premieres et dernieres colonnes N_pts_plus fois.

			for (int j = 0; j < N_pts_plus - 1; j++)
			{
				dec = N0*j;
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i + dec] = M_t_t0[i];
				}

			}

			for (int j = 0; j < N; j++)
			{
				dec = N0 * (N_pts_plus - 1 + j);
				dec1 = N0*j;

				for (int i = 0; i < N0; i++)
				
				{
					M_t_t0_loc[i + dec] = M_t_t0[i + dec1];
				}
			}

			dec1 = N0 * (N - 1);
			for (int j = 0; j < N_pts_plus; j++)
			{
				dec = N0 * (N_pts_plus + N - 1 + j);
				
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i + dec] = M_t_t0[i +dec1 ];
				}

			}
			//A partir de M_t_t0_loc, additionner les matrices d�cal�es :

			//memset a zero car M_t_t0 et M_decalage_t0 sont les memes pointeurs que je donne en input
			//a enlever si nuisible 

			memset(M_decalage_t0, 0, sizeof(double) * N0 * N);

			for (int indice_de_N_dt0_s_dt = 0; indice_de_N_dt0_s_dt < N_dt0_s_dt; indice_de_N_dt0_s_dt++)
			{
				for (int j = 0; j < N; j++)
				{
					dec = N0 * j;
					dec1 = N0 *(indice_de_N_dt0_s_dt + j);
					for (int i = 0; i < N0; i++)
					{
						M_decalage_t0[i + dec] += M_t_t0_loc[i + dec1];
					}
				}
			}
		}

		if (!(N_dt0_s_dt % 2))
		{
			M_t_t0_loc = malloc(sizeof(double)*N0*(N + 2 * N_pts_plus));
			//cr�er d'abord M_t_t0_loc en ajoutant les premieres et dernieres colonnes N_pts_plus fois.

			//ajout de N colonnes � gauche et N colonnes � droite de la matrice donn�e en input
			for (int j = 0; j < N_pts_plus; j++)
			{
				dec = N0 * j;
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i + dec] = M_t_t0[i];
				}

			}

			for (int j = 0; j < N; j++)
			{
				dec = N0 * (N_pts_plus + j);
				dec1 = N0 * j;
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i +dec ] = M_t_t0[i +dec1 ];
				}
			}

			dec1 = N0 * (N - 1);
			for (int j = 0; j < N_pts_plus; j++)
			{
				dec = N0 * (N_pts_plus + N + j);

				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i + dec] = M_t_t0[i + dec1 ];
				}

			}
			//A partir de M_t_t0_loc, additionner les matrices d�cal�es :

			//memset a zero car M_t_t0 et M_decalage_t0 sont les memes pointeurs que je donne en input
			//a enlever si nuisible 
			memset(M_decalage_t0, 0, sizeof(double) * N0 * N);

			for (int indice_de_N_dt0_s_dt = 0; indice_de_N_dt0_s_dt < N_dt0_s_dt; indice_de_N_dt0_s_dt++)
			{
				for (int j = 0; j < N; j++)
				{
					dec = N0 * j;
					dec1 = N0 * (indice_de_N_dt0_s_dt + j);
					for (int i = 0; i < N0; i++)
					{
						M_decalage_t0[i + dec] += + M_t_t0_loc[i + dec1];
					}
				}
			}

		}
		free(M_t_t0_loc);
	}

	if (IM_M_t_t0 != NULL)
	{
		double * M_t_t0_loc;
		double * IM_M_t_t0_loc;
		if (N_dt0_s_dt % 2)
		{
			M_t_t0_loc = malloc(sizeof(double)*N0*(N + 2 * N_pts_plus - 1));
			IM_M_t_t0_loc = malloc(sizeof(double)*N0*(N + 2 * N_pts_plus - 1));

			//cr�er d'abord M_t_t0_loc en ajoutant les premieres et dernieres colonnes N_pts_plus fois.

			for (int j = 0; j < N_pts_plus - 1; j++)
			{
				dec = N0 * j;
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i +dec ] = M_t_t0[i];
					IM_M_t_t0_loc[i + dec] = IM_M_t_t0[i];
				}

			}

			for (int j = 0; j < N; j++)
			{
				dec = N0 * (N_pts_plus - 1 + j);
				dec1 = N0 * j;
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i +dec ] = M_t_t0[i +dec1 ];
					IM_M_t_t0_loc[i + dec] = IM_M_t_t0[i + dec1];
				}
			}

			dec1 = N0 * (N - 1);
			for (int j = 0; j < N_pts_plus; j++)
			{
				dec = N0 * (N_pts_plus + N - 1 + j);
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i +dec] = M_t_t0[i + dec1];
					IM_M_t_t0_loc[i + dec] = IM_M_t_t0[i + dec1];

				}

			}

			//A partir de M_t_t0_loc, additionner les matrices d�cal�es :

			memset(M_decalage_t0, 0, sizeof(double) * N0 * N);
			memset(IM_M_decalage_t0, 0, sizeof(double) * N0 * N);

			for (int indice_de_N_dt0_s_dt = 0; indice_de_N_dt0_s_dt < N_dt0_s_dt; indice_de_N_dt0_s_dt++)
			{
				for (int j = 0; j < N; j++)
				{
					dec = N0 * j;
					dec1 = N0 * indice_de_N_dt0_s_dt + N0 * j;
					for (int i = 0; i < N0; i++)
					{
						M_decalage_t0[i +dec ] +=  M_t_t0_loc[i + dec1];
						IM_M_decalage_t0[i + dec] += IM_M_t_t0_loc[i + dec1];

					}
				}
			}
		}

		if (!(N_dt0_s_dt % 2))
		{
			M_t_t0_loc = malloc(sizeof(double)*N0*(N + 2 * N_pts_plus));
			IM_M_t_t0_loc = malloc(sizeof(double)*N0*(N + 2 * N_pts_plus));

			//cr�er d'abord M_t_t0_loc en ajoutant les premieres et dernieres colonnes N_pts_plus fois.

			for (int j = 0; j < N_pts_plus; j++)
			{
				dec = N0 * j;
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i +dec] = M_t_t0[i];
					IM_M_t_t0_loc[i + dec] = IM_M_t_t0[i];

				}

			}

			for (int j = 0; j < N; j++)
			{
				dec = N0 * (N_pts_plus + j);
				dec1 = N0 * j;
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i +dec ] = M_t_t0[i + dec1];
					IM_M_t_t0_loc[i + dec] = IM_M_t_t0[i + dec1];

				}
			}
			
			dec1 = N0 * (N - 1);
			for (int j = 0; j < N_pts_plus; j++)
			{
				dec = N0 * (N_pts_plus + N + j);
				for (int i = 0; i < N0; i++)
				{
					M_t_t0_loc[i +dec ] = M_t_t0[i + dec1];
					IM_M_t_t0_loc[i +dec] = IM_M_t_t0[i + dec1];

				}

			}
			
			//A partir de M_t_t0_loc, additionner les matrices d�cal�es :

			memset(M_decalage_t0, 0, sizeof(double) * N0 * N);
			memset(IM_M_decalage_t0, 0, sizeof(double) * N0 * N);

			for (int indice_de_N_dt0_s_dt = 0; indice_de_N_dt0_s_dt < N_dt0_s_dt; indice_de_N_dt0_s_dt++)
			{
				for (int j = 0; j < N; j++)
				{
					 dec= N0 * j;
					 dec1 = N0 * indice_de_N_dt0_s_dt + dec;
					for (int i = 0; i < N0; i++)
					{
						M_decalage_t0[i +dec ] +=  M_t_t0_loc[i +dec1 ];
						IM_M_decalage_t0[i + dec] += IM_M_t_t0_loc[i + dec1];

					}
				}
			}

		}
		free(M_t_t0_loc);
		free(IM_M_t_t0_loc);
	}

		
	return;
}

/*void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs != 1 || nrhs != 3)
	{
		mexErrMsgTxt("Wrong number of arguments during call");
		return;
	}
	
	double* IM_M_t_t0;
	if (mxIsComplex(prhs[1]))
	{
		IM_M_t_t0 = mxGetPi(prhs[1]);
	}
	
	else
	{
		IM_M_t_t0 = NULL;
	}


	int m = mxGetM(prhs[0]); //1;
	if (m != 1)
	{
		mexErrMsgTxt("Wrong type of input, first entry must me be a row vector");
		return;
	}
	int mm= mxGetM(prhs[2]);
	if (mm != 1)
	{
		mexErrMsgTxt("Wrong type of input, third entry must me be a row vector");
		return;
	}
	int N0 = mxGetN(prhs[0]);//207
	int N = mxGetN(prhs[2]); //201;


	//INPUTS
	double*t0 = mxGetPr(prhs[0]);
	double *M_t_t0 = mxGetPr(prhs[1]);	
	double*t = mxGetPr(prhs[2]);
	double dt, dt0;
	dt = t[9] - t[8]; dt0 = t0[9] - t0[8];
	int N_dt0_s_dt = dt0 / dt;
	int N_pts_plus = (int)(floor((N_dt0_s_dt + 1.0) / 2.0));

	
	//OUTPUTS
	double*M_decalage_t0;
	double*IM_M_decalage_t0;
	if (IM_M_t_t0 == NULL)
	{
		IM_M_decalage_t0 = NULL;
		plhs[0] = mxCreateDoubleMatrix(N0, N, mxREAL);
		M_decalage_t0 = mxGetPr(plhs[0]);
		if (N_dt0_s_dt == 0)
		{
			for (int i = 0; i < N0; i++)
			{
				for (int j = 0; j < N; j++)
				{
					M_decalage_t0[i + N0 * j] = M_t_t0[i + N0 * j];
				}
			}

			return;
		}
	}
	if (IM_M_t_t0 != NULL)
	{
		plhs[0] = mxCreateDoubleMatrix(N0, N, mxCOMPLEX);
		M_decalage_t0 = mxGetPr(plhs[0]);
		IM_M_decalage_t0 = mxGetPi(plhs[0]);
		if (N_dt0_s_dt == 0)
		{
			for (int i = 0; i < N0; i++)
			{
				for (int j = 0; j < N; j++)
				{
					M_decalage_t0[i + N0 * j] = M_t_t0[i + N0 * j];
					IM_M_decalage_t0[i + N0 * j] = IM_M_t_t0[i + N0 * j];
				}
			}

			return;
		}
	}
	
	
	fc_somme_sur_les_dt0_decales_Matrice_de_t_decalee_de_t0_BEN(N_pts_plus,N_dt0_s_dt, M_t_t0, IM_M_t_t0, M_decalage_t0, IM_M_decalage_t0, 
		N0,N);
	return;
}*/