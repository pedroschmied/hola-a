#include "interaccion.h"
#include "general.h"
#include "avanzar.h"

double velocity_verlet(double *x, double *v, double *f, int N, double h, double L, double *tabla_F, double *tabla_V, double dr2)
{	
	double rc2 = 6.25, Vij;
	position_verlet(x, v, N, h, f);
	apply_PBC(x, N, L);
	velocidad(v, N, h, f);	
	double Vij = fuerzas(tabla_F, tabla_V, f, x, rc2, dr2, N, L);
	velocidad(v, N, h, f);
	return Vij;
}
int apply_PBC(double *x, int N, double L)
{
	int i;
	for(i = 0; i < 3 * N; i ++)
	{
		if (*(x + i) > L)
		{
			*(x +i) = *(x +i) - L;
		}
		else if (*(x + i) < 0)
		{
			*(x +i) = *(x +i) + L;
		}
	}
	return 0;
}

// h: paso temporal
int position_verlet(double *x, double *v, int N, double h, double *F)
{
	int i, k;
	for(i = 0; i < N; i ++)
	{
		for(k = 0; k < 3; k++)
		{
			*(x + 3 * i + k) += *(v + 3 * i + k) * h + *(F + 3 * i + k) * h * h / 2.0; //x(t + h) = x(t) + v(t)*h + f(t)/(2*m) * h^2;
		}
	}
	return 0;
}
int velocity_verlet(double *v, int N, double h, double *F)
{
	int i, k;
	for(i = 0; i < N; i ++)
	{
		for(k = 0; k < 3; k++)
		{
			*(v + 3 * i + k) +=  (*(F + 3 * i + k) * h / 2.0; 
		}
	}
	return 0;
}
int fuerza_PCB(double *delta_r, double L)
{
	int k;
	for (k = 0; k < 3; k++)
	{
		if(*(delta_r + k) > L / 2.0)
		{
			*(delta_r + k) -= L;
		}
		else if(*(delta_r + k) <= - L / 2.0)
		{
			*(delta_r + k) += L;
		}
	}
	return 0;
}
double fuerzas(double *tabla_F, double *tabla_V, double *F, double *x, double rc2, double dr2, int N, double L)
{
	double  *F_mod;
	F_mod = (double*) malloc(1 * sizeof(double));
	double *delta_r;
	delta_r = (double*) malloc(3 * sizeof(double));
	double A, rij2, Vij, potencial = 0.0;
	int i, j, k;
	for (i = 0; i < 3 * N; i++)
	{
		*(F2 + i) = *(F + i);
		*(F + i) = 0.0;
	}
	for (i = 0; i < N - 1; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			delta_x(x, i, j, delta_r);
			fuerza_PCB(delta_r, (float)L);
			rij2 = norma2(delta_r);
			if(rij2 < rc2)
			{
				Vij = pair_force(tabla_F, tabla_V, rij2, dr2, F_mod);
//				r = sqrt(rij2);
				for (k = 0; k < 3; k++)
				{
					A = *F_mod * *(delta_r + k);
					*(F + 3 * i + k) += A;
					*(F + 3 * j + k) -= A;
				}
				potencial += Vij;
			}
		}
	}
	free(F_mod);
	free(delta_r);
	return potencial;
}
