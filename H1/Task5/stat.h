
#ifndef _stat_h
#define _stat_h

double get_mean(double *vector, int sz);
double get_variance(double *vector, double mu, int sz);

double get_mean_norm(double (*vector)[3], int sz);
double get_variance_norm(double (*vector)[3], double mean, int sz);

double velocityCorr(double (*vel)[3], double (*prevVel)[3], int N);
double msdCorr(double (*pos)[3], double (*prevPos)[3], int N) ;

#endif