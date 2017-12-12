
#ifndef _h2_h
#define _h2_h

double energy2(double R1[3], double R2[3], double alpha);
double energy(double R[6], double alpha);

double absWaveFunction(double * R, double alpha);
double trialWaveFunction(double * R, double alpha);
double probFunction(double * Rnew, double * Rold, double alpha);

#endif
