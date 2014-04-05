/* bubble.h
 *
 * Code relating to the calculation of the bubble power spectrum
 * My crude halo based model
 */

#ifndef BUBBLE_H
#define BUBBLE_H

double uTOP(double k, double r);
double PS1bxx1V(double k, double z, CosmParam *c);
double PS1bxx(double k, double z, CosmParam *c);
void setPS1bxxIntegrand(double m, double y[], double deriv[], double k1, double z1, Cosmology *c1, int flag);
void ps1bxxIntegrand(double m, double y[], double deriv[]);
//filling fraction
double fillFrac(double z, CosmParam *c);
void setFillFracIntegrand(double m, double y[], double deriv[], double k1, double z1, Cosmology *c1, int flag);
void fillFracIntegrand(double m, double y[], double deriv[]);

#endif
