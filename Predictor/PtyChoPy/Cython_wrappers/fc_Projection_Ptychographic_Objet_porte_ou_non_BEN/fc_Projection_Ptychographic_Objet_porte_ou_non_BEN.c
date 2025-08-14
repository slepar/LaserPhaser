/*********************************************************************************************************
* Date: 6 juillet 2018																					 *
* Auteur: Benoit Brizard																				 *
*																										 *
* La projection ptychographic calcule le produit temporel des fonctions									 *
* decalees � l'iteration i+2 en fonction de celui-ci � literation i+1 en								 *
* appliquant la contrainte ptychographique :															 *
* Apr�s retour dans l'espace temporel (apres la contrainte XP), on isole								 *
* les 2 fonctions objet et champ par la contrainte ptycho, pui on calcule							     *
* de nouveau ce m�me produit temporel, mais uniquement avec ces 2 fonctions								 *
* isol�e.																							     *
* De plus, la fc renvoie les deux fonction temporelles  objet et champ									 *
* R�f�rence : page 147 th�se d'Adrien Leblanc, correspondance avec la variable Pi_OS					 *
* INPUTS : Produit_t_t0,  t0,  t, k_max  , Case_objet													 *
* OUTPUTS : M_Produit_t_t0_it_plus2 , Obj_t, Chp_t														 *
*********************************************************************************************************/

#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
// #include <mex.h>
#include "fc_decalage_temporelle_vecteur_t_de_t0_BEN.h"
#include "fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN.h"
#include "fc_Projection_Ptychographic_Objet_porte_ou_non_BEN.h"

//On a besoin dans ce programme de la fonction "circshift", qui est une fonction "built-in" de Matlab
//Ci-bas on retrouve l'algorithme original, puis une version l�g�rement modifi�e pour l'usage sp�cifique qu'on
//en fait ici


/*
*
* Modulo are NEVER negative, so it is the trick needed for circshift
* voici l'impl�mentation d'un modulo jj = mod(a,b) : 
* int jj = a % b;
* if (jj<0) jj = b + jj;
* retur jj;
*/

// VOICI l'algorithme circshift ORIGINAL, NOUS ON SHIFT RIEN EN Y, d'o� ma petite modification de la fonction (voir plus bas, apr�s la section mise en commentaire)
/*void circshift(double * REAL_out, double * IM_out, double * REAL_in, double * IM_in, int xdim, int ydim, int xshift, int yshift) 
{
	for (int i = 0; i < xdim; i++) 
	{
		int ii = (i + xshift) % xdim;
		if (ii<0) ii = xdim + ii;
		for (int j = 0; j < ydim; j++) 
		{
			int jj = (j + yshift) % ydim;
			if (jj<0) jj = ydim + jj;
			REAL_out[ii * ydim + jj] = REAL_in[i * ydim + j];
			IM_out[ii * ydim + jj] = IM_in[i * ydim + j];
		}
	}
	return;
}*/

void circshift(double * REAL_out, double * IM_out, double * REAL_in, double * IM_in, int xdim, int ydim, int xshift) //xdim=N, ydim=N0
{
	for (int i = 0; i < xdim; i++)
	{
		int ii = (i + xshift) % xdim;
		if (ii<0) ii = xdim + ii;
		for (int j = 0; j < ydim; j++)
		{
			REAL_out[ii * ydim + j] = REAL_in[i * ydim + j];
			IM_out[ii * ydim + j] = IM_in[i * ydim + j];
		}
	}
	return;
}

void fc_Projection_Ptychographic_Objet_porte_ou_non_BEN(double max_t0, double dt0, double dt, double * t, int N_pts_plus_somme, int N_pts_plus_decalage, int N_dt0_s_dt,
	double k_max, double case_object, double *REAL_Produit_t_t0, double *IM_Produit_t_t0, int N0, int N, double *REAL_Obj_t, double *IM_Obj_t, double *REAL_Chp_t,
	double * IM_Chp_t, double*REAL_M_Produit_t_t0_it_plus2, double*IM_M_Produit_t_t0_it_plus2)
{
	// Contrainte ptychographic :
	// On deduit les functions Objet et Champ qui correspondent au mieux le produit d�cal�

	fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(max_t0, dt0, dt, t, N_pts_plus_somme, N_pts_plus_decalage, N_dt0_s_dt, k_max, case_object,
		REAL_Produit_t_t0, IM_Produit_t_t0, N0, N, REAL_Obj_t, IM_Obj_t, REAL_Chp_t, IM_Chp_t);
	
	double sign = 1.0, max = 0.0,current=0.0;
	int dec, somme, indice;

	// Nouveau produit d�cal� construit � partir des 2 fonctions retrouv�es
	// on decale l'objet de chaque t0
	fc_decalage_temporelle_vecteur_t_de_t0_BEN(dt0, dt, N_pts_plus_decalage, REAL_Obj_t, IM_Obj_t, &sign, N0, N, REAL_M_Produit_t_t0_it_plus2, IM_M_Produit_t_t0_it_plus2);
	
	// Maintenant on effectue produit temporel du champ et de l'objet d�cal�.
	// Cela impliquera la multiplication de nombres complexes.

	// ici on va faire la multiplication de 2 nombres complexes,  on a de fa�on g�n�rale :
	// (a+i*b)*(c+d). Par souci d'optimisation, on peut trouver le r�sultat avec seulement 3 multiplications. 
	// La partie r�elle du totale est : a*c-b*d. La partie imaginaire du totale est : (a+c)*(c+d) -a*c-b*d ->la variable contenant ce r�sultat est nomm�e "big".
	// C'est ce qu'on impl�mente ici.

	if (case_object != 1)//on multiplie 2 complex
	{
		double ac = 0, bd = 0, big = 0;
		for (int j = 0; j < N; j++)
		{
			dec = N0*j;
			for (int i = 0; i < N0; i++)
			{
				somme = dec + i;
				ac = REAL_M_Produit_t_t0_it_plus2[somme] * REAL_Chp_t[j];
				bd = IM_M_Produit_t_t0_it_plus2[somme] * IM_Chp_t[j];
				big = (REAL_M_Produit_t_t0_it_plus2[somme]+ IM_M_Produit_t_t0_it_plus2[somme])*(REAL_Chp_t[j] + IM_Chp_t[j]) - ac - bd;
				REAL_M_Produit_t_t0_it_plus2[somme] = ac - bd;
				IM_M_Produit_t_t0_it_plus2[somme] = big;

			}
		}

	}
	if (case_object == 1)//multiplie un reel et un complex 
	{
		for (int j = 0; j < N; j++)
		{
			dec = N0*j;
			for (int i = 0; i < N0; i++)
			{
				somme = dec + i;
				IM_M_Produit_t_t0_it_plus2[somme] = IM_Chp_t[j] * REAL_M_Produit_t_t0_it_plus2[somme];
				REAL_M_Produit_t_t0_it_plus2[somme] *= REAL_Chp_t[j];


			}
		}
		
	}
	

	// Enlever grossi�rement la phase lin�aire en recentrant temporellement le produit � t = 0
	
	for (int j = 0; j < N;j++)
	{
		dec = j*N0;
		for (int i = 0; i < N0; i++)
		{
			somme = dec + i;
			current += REAL_M_Produit_t_t0_it_plus2[somme] * REAL_M_Produit_t_t0_it_plus2[somme] + IM_M_Produit_t_t0_it_plus2[somme] * IM_M_Produit_t_t0_it_plus2[somme];
		}
		if (current > max)
		{
			max = current;
			indice = j;
		}
		current = 0.0;
	}
	
	/*
	* En g�n�ral, on a jamais besoin de faire le circshift, soit  (floor(N / 2) - indice - 1)=0. Afin d'optimiser le code,
	* jamais on call circshift, sauf lorsque n�cessaire. C'est alors qu'on alloue de la m�moire. C'est plus court overall 
	* comme �a
	*/
	

	indice = (int)(floor(N / 2) - (indice));
	if (dec != 0)
	{
		double*REAL_temp = malloc(sizeof(double)*N*N0);
		double*IM_temp = malloc(sizeof(double)*N*N0);
		for (int i = 0; i < N; i++)
		{
			dec = N0*i;
			for (int j = 0; j < N0; j++)
			{
				somme = dec + j;
				REAL_temp[somme] = REAL_M_Produit_t_t0_it_plus2[somme];
				IM_temp[somme] = IM_M_Produit_t_t0_it_plus2[somme];
			}
		}
		circshift(REAL_M_Produit_t_t0_it_plus2, IM_M_Produit_t_t0_it_plus2, REAL_temp, IM_temp, N, N0, indice);
		free(REAL_temp);
		free(IM_temp);
	}

	return;

}





// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// {
// 	if (nlhs != 3 || nrhs != 5)
// 	{
// 		mexErrMsgTxt("Wrong number of arguments during call");
// 	}
	
// 	if (!mxIsComplex(prhs[0]))
// 	{
// 		mexWarnMsgTxt("This program assumes that first input is complex.");
// 	}
	
// 	int N0 = mxGetM(prhs[0]); //207;
// 	int N = mxGetN(prhs[0]); //201;


// 							 /*
// 							 double *Ar = mxGetPr(A IN); // Real data
// 							 double *Ai = mxGetPi(A IN); // Imaginary data
// 							 And Ar[m + M * n] and Ai[m + M * n] are the real and imaginary parts of A(m + 1, n + 1).
// 							 To create a 2 - D complex array, use
// 							 B OUT = mxCreateDoubleMatrix(M, N, mxCOMPLEX);*/

// 							 //INPUTS
// 	double *REAL_Produit_t_t0 = mxGetPr(prhs[0]);
// 	double *IM_Produit_t_t0 = mxGetPi(prhs[0]);
// 	double*t0 = mxGetPr(prhs[1]);
// 	double*t = mxGetPr(prhs[2]);
// 	double*k_max = mxGetPr(prhs[3]);
// 	double*case_object = mxGetPr(prhs[4]);
// 	double dt, dt0, max_t0;
// 	max_t0 = t0[0];
// 	for (int i = 1; i < N0; i++)
// 	{
// 		if (t0[i] > max_t0)
// 			max_t0 = t0[i];
// 	}
// 	dt = t[9] - t[8]; dt0 = t0[9] - t0[8];
// 	int N_dt0_s_dt = dt0 / dt;
// 	int N_pts_plus_somme = (int)(floor((N_dt0_s_dt + 1.0) / 2.0));//somme decale
// 	int N_pts_plus_decalage = (int)((N0 - 1)* N_dt0_s_dt);//decalge temporelle MATRICE et VECTEUR


// 		//OUTPUTS
// 	plhs[0] = mxCreateDoubleMatrix(N0, N, mxCOMPLEX);
// 	double*REAL_M_Produit_t_t0_it_plus2= mxGetPr(plhs[0]);
// 	double*IM_M_Produit_t_t0_it_plus2 = mxGetPi(plhs[0]);
// 	plhs[1] = mxCreateDoubleMatrix(1, N, mxCOMPLEX);
// 	double *REAL_Obj_t = mxGetPr(plhs[1]);
// 	double *IM_Obj_t = mxGetPi(plhs[1]);
// 	plhs[2] = mxCreateDoubleMatrix(1, N, mxCOMPLEX);
// 	double *REAL_Chp_t = mxGetPr(plhs[2]);
// 	double *IM_Chp_t = mxGetPi(plhs[2]);
	
// 	fc_Projection_Ptychographic_Objet_porte_ou_non_BEN(max_t0, dt0, dt, t, N_pts_plus_somme, N_pts_plus_decalage, N_dt0_s_dt, *k_max, *case_object,
// 		REAL_Produit_t_t0, IM_Produit_t_t0, N0, N, REAL_Obj_t, IM_Obj_t, REAL_Chp_t, IM_Chp_t, REAL_M_Produit_t_t0_it_plus2, IM_M_Produit_t_t0_it_plus2);
// 	return;
// }