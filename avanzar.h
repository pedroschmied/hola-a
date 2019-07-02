#ifndef AVANZAR_H
#define AVANZAR_H
#include <math.h>
double velocity_verlet(double *x, double *v, double *f, int N, double h, double L, double *tabla_F, double *tabla_V, double dr2);
int apply_PBC(double *x, int N, double L);
int position_verlet(double *x, double *v, int N, double h, double *F);
int velocidad(double *v, int N, double h, double *F);
int fuerza_PCB(double *delta_r, double L);
double fuerzas(double *tabla_F, double *tabla_V, double *F, double *x, double rc2, double dr2, int N, double L);
#endif
