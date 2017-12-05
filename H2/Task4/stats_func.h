
#ifndef _stats_func_h
#define _stats_func_h

extern double getMean(double *, int n);
extern double getSquaredValueMean(double * values, int nbr_of_values);
extern double getAheadMean(double * values, int nbr_of_values, int k);
extern double autoCorrelation(double * values, int nbr_of_values, int k);
extern double d_abs(double d);
extern int findS(double * values, int nbr_of_values, double threshold);
extern double blockAverageS(double * values, int nbr_of_values, int B);
extern double blockAverage(double * values, int nbr_of_values, int j, int B);
extern double sampleVariance(double * values, int nbr_of_values);
#endif
