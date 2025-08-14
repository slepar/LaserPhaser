/*********************************************************************************************************
* Date: 6 juillet 2018																					 *
* Auteur: Benoit Brizard																				 *
*																										 *
* Fonction appliquant la contrainte ptycographique � la matrice Produit_t_t0. Cette fonction est		 *
* souvent appel�e par la fonction plus g�n�rale fc_Projection_Ptychographic_Objet_porte_ou_non_BEN		 *
* Pour les d�tails, se r�f�rer au script Matlab et � la th�se d'Adrien Leblanc. � la page 230, on a	les	 *
* �quivalences suivantes pour passer de la notation de la th�se au nom des variables de ce programme :   *
* a_i-1(t,t0) = Produit_t_t0																			 *
* a_i-1(t+t0,t0) = Produit_t_t0_decal																	 *
* INPUTS : t0, Produit_t_t0, t, k_max , Case_objet													     *
* OUTPUTS : Obj_t, Chp_t																				 *
*********************************************************************************************************/


#define _CRT_SECURE_NO_WARNINGS
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
// #include <mex.h>
#include "fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN.h"
#include "fc_decalage_temporelle_Matrice_t_de_t0_BEN.h"
#include "fc_decalage_temporelle_vecteur_t_de_t0_BEN.h"
#include "fc_somme_sur_les_dt0_decales_Matrice_de_t_decalee_de_t0_BEN.h"

										/*********************************************************************************************************
										*	   SOUS-FONCTIONS N�CESSAIRES POUR FAIRE LE PONT ENTRE DES COMMANDES MATLAB ET LE LANGUAGE EN C :	 *
										*********************************************************************************************************/




/*********************************************************************************************************
* sert � faire l'�quivalent de (vient de Matlab) :														 *
*	 Chp_t_k = sum(conj(M_Obj_decalage_Sum).*Produit_t_t0, 1)  . / sum(abs(M_Obj_decalage_Sum). ^ 2, 1); *
*	 Chp_t_k(isnan(Chp_t_k)) = 0;																		 *
* en une fonction en C																					 *
* INPUTS : M_Obj_decalage_Su,, IM_M_Obj_decalage_SUM,REAL_Produit_t_t0_decal, IM_Produit_t_t0_decal et N0*
* et N,	qui sont les dimensions matricielles, respectivement le nombre de lignes et le nombre de colonnes*
* OUTPUTS :  REAL_Obj_t, IM_Obj_t																		 *
*********************************************************************************************************/
void sum_colonne_normale(double*REAL_Obj_t, double*IM_Obj_t, double*M_Obj_decalage, double*IM_M_Obj_decalage, double*REAL_Produit_t_t0_decal,
	double*IM_Produit_t_t0_decal,int N0, int N)
{	
	double * diviseur = calloc(N,sizeof(double));
	//typiquement 4 cas : 
	//Obj_decalage complexe, Produit_t_t0_decal reel//Obj_decalage complexe, Produit_t_t0_decal complexe
	//Obj_decalage reel,Produit_t_t0_decal reel//Obj_decalage reel, Produit_t_t0_decal complexe

	//les 1er et 4e cas sont impossibles (ne peuvent jamais survenir)

	int indice = 0;
	// 2e cas : soit Obj_decalage complexe, Produit_t_t0_decal complexe
	// ici on va faire la multiplication de 2 nombres complexes, puisqu'on prend le conjugu� du premier nombre, on a de fa�on g�n�rale :
	// (a-i*b)*(c+d). Par souci d'optimisation, on peut trouver le r�sultat avec seulement 3 multiplications. 
	// La partie r�elle du totale est : a*c+b*d. La partie imaginaire du totale est : (a-c)*(c+d) -a*c+b*d ->la variable contenant ce r�sultat est nomm�e "big".
	// C'est ce qu'on impl�mente ici.
	
	if ( IM_Produit_t_t0_decal != NULL)
	{
		double ac = 0, bd = 0, big = 0;
		for (int j = 0; j < N;j++)
		{
			indice = N0 * j;
			REAL_Obj_t[j] = 0.0; IM_Obj_t[j] = 0.0; //on  reset l'array, on repart a neuf (contenait initialement des donn�es)
			for (int i = 0; i < N0;i++)
			{
				ac = M_Obj_decalage[i + indice] * REAL_Produit_t_t0_decal[i + indice];
				bd = IM_M_Obj_decalage[i + indice] * IM_Produit_t_t0_decal[i + indice];
				big = (M_Obj_decalage[i + indice] - IM_M_Obj_decalage[i + indice])*(REAL_Produit_t_t0_decal[i + indice] + IM_Produit_t_t0_decal[i + indice]) - ac + bd;
				REAL_Obj_t[j] = REAL_Obj_t[j] +ac+bd;
				IM_Obj_t[j] = IM_Obj_t[j] + big;
				diviseur[j] = diviseur[j] + M_Obj_decalage[i + indice]* M_Obj_decalage[i + indice] + IM_M_Obj_decalage[i + indice]* IM_M_Obj_decalage[i + indice];
				
			}
					
			REAL_Obj_t[j] /=  diviseur[j]; IM_Obj_t[j] /= diviseur[j];
			if (isnan(REAL_Obj_t[j]) || isnan(IM_Obj_t[j]))
			{
				REAL_Obj_t[j] = 0;
				IM_Obj_t[j] = 0;
			}
		}
		free(diviseur);
		return;
	}

	//3e cas : les 2 matrices � multipier sont r�elles
	if ( IM_Produit_t_t0_decal == NULL)
	{
		for (int j = 0; j < N;j++)
		{
			indice = N0*j;
			REAL_Obj_t[j] = 0.0; //on  reset l'array, on repart a neuf (contenait initialement des donn�es)
			for (int i = 0; i < N0;i++)
			{
				REAL_Obj_t[j] = REAL_Obj_t[j] + M_Obj_decalage[i + indice] * REAL_Produit_t_t0_decal[i + indice];
				diviseur[j] = diviseur[j] + M_Obj_decalage[i + indice] * M_Obj_decalage[i + indice];
				
			}
			
			REAL_Obj_t[j] /= diviseur[j]; 
			
		}
		free(diviseur);
		return;
	}

	free(diviseur);
}

/*********************************************************************************************************
* sert � faire l'�quivalent de (vient de Matlab) :														 *
*	 Chp_t_k = Chp_t_k/mean(abs(Chp_t_k)) .* V_filtre;													 *
* en une fonction en C																					 *
* La variable "version" est utilis�e car au fil du code on a besoin de faire une op�ration TR�S similaire*
* d�pendamment de la variable "case_object". On a donc diff�rente "version" de la m�me fonction de base  *
*********************************************************************************************************/

void mean_abs_filtre(double*r_input, double*im_input, double*V_filtr,int N,double version)
{
	double mean = 0;
	double div = 0;
	if (version == 0.0)
	{
		if (im_input != NULL)
		{
			for (int i = 0; i < N; i++)
			{
				mean += sqrt(r_input[i]* r_input[i] + im_input[i]* im_input[i]);
			}
			mean /=  N;
			for (int i = 0; i < N; i++)
			{
				div= V_filtr[i] / mean;
				r_input[i] *= div;
				im_input[i] *= div;
			
			}
			return;
		}
		if (im_input == NULL)
		{
			for (int i = 0; i < N; i++)
			{
				if (r_input[i] < 0)
				{
					mean -= r_input[i];
				}
				else
				{
					mean += r_input[i];
				}
				
			}
			mean /= N;
			for (int i = 0; i < N; i++)
			{
				r_input[i] *= V_filtr[i]/mean;
			}


		}
		return;
	}

	if (version == 1.0)
	{
		if (im_input != NULL)
		{
			for (int i = 0; i < N; i++)
			{
				r_input[i]= sqrt(r_input[i]* r_input[i] + im_input[i]* im_input[i]);
				mean += r_input[i];
			}
			mean /= N;
			for (int i = 0; i < N; i++)
			{
				r_input[i] *= V_filtr[i]/mean;
				im_input[i] = 0.0;
			}
			return;
		}
		if (im_input == NULL)
		{
			for (int i = 0; i < N; i++)
			{
				if (r_input[i] < 0)
				{
					r_input[i] *= -1.0;
				}
				mean += r_input[i];
			}
			mean /= N;
			for (int i = 0; i < N; i++)
			{
				r_input[i] *= V_filtr[i]/mean;
			}
		}
		return;
	}

}

void fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(double max_t0, double dt0, double dt, double * t, int N_pts_plus_somme, int N_pts_plus_decalage, int N_dt0_s_dt,
	double k_max, double case_object, double *REAL_Produit_t_t0, double *IM_Produit_t_t0, int N0, int N, double *REAL_Obj_t,
	double *IM_Obj_t, double *REAL_Chp_t, double * IM_Chp_t)
{
	double * REAL_Produit_t_t0_decal = malloc(sizeof(double)*N0*N);
	double * V_filtre = malloc(sizeof(double)*N);


	double Ordre_filtr = 40.0;//malheureusement pour imiter la fonction parfaitement...ces chiffres sont donn�s directement,
	double v_max_t0 = 0.9999; //ce ne sont pas des inputs de la fonction originale en Matlab. Il sera difficile de changer ces param�tres
							  //une fois que le fichier ex�cutable sera cr��...a penser de cr�er des variables qu'on peut modifier
							  //pour ces parametres (si pertinent, ca peut �tre utile de changer ces facteurs qui m'ont l'air d'avoir
							  //�t� attribu� "a l'oeil" sans fondement theorique.

	/* determiner zone ou on considere la reconstruction : autour du temps du scan
		% je fais un filtre hypergaussien d'ordre ''Ordre_filtr'', tel que sa valeur au
		% max de t0 soit ''v_max_t0'' tres proche de 1, on calcule alors la
		% largeur du filtre adequat ''tau_filtr = t0 / (-ln(v)) ^ (1 / ordre) ''
		Ordre_filtr = 40;
		v_max_t0 = 0.9999;  % valeur du filtre au niveau des bords du scan*/
	double to_pow = max_t0/pow(-1.0*log(v_max_t0), (1.0 / Ordre_filtr));
							 // d�duction de la larguer du filtre
							// to_pow est maintenant l'�quivalent de tau_filtr du code matlab

	//on profite de la meme boucle pour initialiser obj_t_k en meme temps de creer le V_filtre
	//srand(time(NULL));
	srand((unsigned int)time(NULL));

	double mean = 0;
	for (int i = 0; i < N; i++)
	{
		V_filtre[i] = exp(-1.0*pow((t[i]) / to_pow, Ordre_filtr));
		REAL_Obj_t[i] = rand() / ((double)RAND_MAX);
		mean += REAL_Obj_t[i];//pour calculer la moyenne 
	}
	mean /= N;
	for (int i = 0; i < N; i++)
	{
		REAL_Obj_t[i] /=  mean;
	}

	double sign = -1.0;
	double k = 1;
	double*IM_Produit_t_t0_decal;

	if (IM_Produit_t_t0 == NULL)
	{
		IM_Produit_t_t0_decal = NULL;
	}
	if (IM_Produit_t_t0 != NULL)
	{
		IM_Produit_t_t0_decal = malloc(sizeof(double)*N0*N);
	}
	fc_decalage_temporelle_Matrice_t_de_t0_BEN(dt0, dt, N_pts_plus_decalage, REAL_Produit_t_t0, IM_Produit_t_t0, &sign, N0, N,
		REAL_Produit_t_t0_decal, IM_Produit_t_t0_decal);


	//M_Obj_decalage servira vraiment de "conteneur" dans lequel je vais sauvegarder toute matrice interm�diaire au calcul
	double*M_Obj_decalage = calloc(N*N0, sizeof(double));
	double*IM_M_Obj_decalage = calloc(N*N0, sizeof(double));
	double buffer;//utilis� pour case_object=2
	while (k <= k_max)
	{
		k++;
		sign = 1.0;
		fc_decalage_temporelle_vecteur_t_de_t0_BEN(dt0, dt, N_pts_plus_decalage, REAL_Obj_t, IM_Obj_t, &sign, N0, N,
			M_Obj_decalage, IM_M_Obj_decalage);


	
		fc_somme_sur_les_dt0_decales_Matrice_de_t_decalee_de_t0_BEN(N_pts_plus_somme, N_dt0_s_dt, M_Obj_decalage, IM_M_Obj_decalage,
			M_Obj_decalage, IM_M_Obj_decalage, N0, N);
		
		//	M_Obj_decalage est maintenant l'�quivalent de M_Obj_decalage_Sum dans le code Matlab 
		
		sum_colonne_normale(REAL_Chp_t, IM_Chp_t, M_Obj_decalage, IM_M_Obj_decalage, REAL_Produit_t_t0,
			IM_Produit_t_t0, N0, N);
		mean_abs_filtre(REAL_Chp_t, IM_Chp_t, V_filtre, N, 0);
		sign = -1.0;
		fc_decalage_temporelle_vecteur_t_de_t0_BEN(dt0, dt, N_pts_plus_decalage, REAL_Chp_t, IM_Chp_t, &sign, N0, N,
			M_Obj_decalage, IM_M_Obj_decalage);

		//	M_Obj_decalage est maintenant l'�quivalent de M_Chp_decalage dans le code Matlab 

		sum_colonne_normale(REAL_Obj_t, IM_Obj_t, M_Obj_decalage, IM_M_Obj_decalage, REAL_Produit_t_t0_decal,
			IM_Produit_t_t0_decal, N0, N);

		if (case_object == 0.0)
		{
			mean_abs_filtre(REAL_Obj_t, IM_Obj_t, V_filtre, N, case_object);
		}
		else if (case_object == 1.0) // Objet reel(chgt de transmission mais pas de phase)
		{
			mean_abs_filtre(REAL_Obj_t, IM_Obj_t, V_filtre, N, case_object);
		}

		else if (case_object == 2.0)// Objet de phase(amplitude ne varie pas)
		{
			buffer = 0;
			for (int i = 0; i < N; i++)
			{
				buffer = sqrt(REAL_Obj_t[i] * REAL_Obj_t[i] + IM_Obj_t[i] * IM_Obj_t[i]);
				REAL_Obj_t[i] /= buffer;
				IM_Obj_t[i] /= buffer;

				if (isnan(REAL_Obj_t[i]) || isnan(IM_Obj_t[i]))
				{
					REAL_Obj_t[i] = 0;
					IM_Obj_t[i] = 0;
				}

			}

		}

	}
	
				
	free(REAL_Produit_t_t0_decal);
	free(V_filtre);
	if (IM_Produit_t_t0_decal != NULL)
	{
		free(IM_Produit_t_t0_decal);
	}
	free(M_Obj_decalage);
	free(IM_M_Obj_decalage);
	return;
}



/*void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs != 2 || nrhs != 5)
	{
		mexErrMsgTxt("Wrong number of arguments during call");
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
	
	//INPUTS

	/*double *REAL_Produit_t_t0 = mxGetPr(prhs[1]);
	double*t0 = mxGetPr(prhs[0]);
	double*t = mxGetPr(prhs[2]);
	double*k_max = mxGetPr(prhs[3]);
	double*case_object = mxGetPr(prhs[4]);
	double dt, dt0, max_t0;
	max_t0 = t0[0];
	for (int i = 1; i < N0; i++)
	{
		if (t0[i] > max_t0)
			max_t0 = t0[i];
	}
	dt = t[9] - t[8]; dt0 = t0[9] - t0[8];
	int N_dt0_s_dt = dt0 / dt;
	int N_pts_plus_somme = (int)(floor((N_dt0_s_dt + 1.0) / 2.0));//somme decale
	int N_pts_plus_decalage = (int)((N0 - 1)* N_dt0_s_dt);//decalge temporelle MATRICE et VECTEUR


	//OUTPUTS

	double *REAL_Obj_t;
	double *IM_Obj_t;
	double *REAL_Chp_t;
	double *IM_Chp_t;
	if (IM_Produit_t_t0 != NULL)
	{
		plhs[0] = mxCreateDoubleMatrix(1, N, mxCOMPLEX);
		REAL_Obj_t = mxGetPr(plhs[0]);
		IM_Obj_t = mxGetPi(plhs[0]);
		plhs[1] = mxCreateDoubleMatrix(1, N, mxCOMPLEX);
		REAL_Chp_t = mxGetPr(plhs[1]);
		IM_Chp_t = mxGetPi(plhs[1]);
	}
	else
	{
		plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
		REAL_Obj_t = mxGetPr(plhs[0]);
		IM_Obj_t = NULL;
		plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
		REAL_Chp_t = mxGetPr(plhs[1]);
		IM_Chp_t = NULL;
	}
	fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(max_t0,dt0, dt,t, N_pts_plus_somme, N_pts_plus_decalage, N_dt0_s_dt,*k_max, *case_object,
		REAL_Produit_t_t0, IM_Produit_t_t0,	N0, N, REAL_Obj_t, IM_Obj_t, REAL_Chp_t, IM_Chp_t);
	return;
}*/