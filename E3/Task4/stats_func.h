
#ifndef _stats_func_h
#define _stats_func_h

extern double getMean(double *, int n);
extern double getSquaredValueMean(double * values, int nbr_of_values);
extern double getAheadMean(double * values, int nbr_of_values, int k);
double autoCorrelation(double * values, int nbr_of_values, int k);
#endif
