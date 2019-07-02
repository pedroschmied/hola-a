#ifndef INICIALIZAR_H
#define INICIALIZAR_H
#include <math.h>
double set_pos(double *x, int N, double L);
double set_vel(double *v, int N, double T);
double inicializar(double *x, double *v, double *f, int N, double L, double T);
#endif
