#include "stats_func.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double getMean(double * values, int nbr_of_values) {
    int i;
    double sum = 0;
    for(i = 0; i < nbr_of_values; i++)
        sum += values[i];
    return sum/nbr_of_values;    
}

double getSquaredValueMean(double * values, int nbr_of_values) {
    int i;
    double sum = 0;
    for(i = 0; i < nbr_of_values; i++)
        sum += values[i]*values[i];
    return sum/nbr_of_values;    
}

double getAheadMean(double * values, int nbr_of_values, int k) {
    if(k < nbr_of_values) {    
        int i;
        double sum = 0;
        for(i = 0; i < nbr_of_values - k; i++) {
            sum += values[i]*values[i+k];
        }
        return sum/(nbr_of_values - k);
    } 
    else { 
           printf("Error: getAheadMean k >= N \n");     
           return 0;
    }
}

double autoCorrelation(double * values, int nbr_of_values, int k) {
    double phi_i = getMean(values, nbr_of_values);
    phi_i *= phi_i;
    double phi_ik = getAheadMean(values, nbr_of_values, k);
    double phiSq = getSquaredValueMean(values, nbr_of_values);
    return (phi_ik - phi_i)/(phiSq - phi_i);
}

int findS(double * values, int nbr_of_values, double threshold) {
    int k;
    int s = -1;     
    double target = 0.135;
    for(k = 0; k < nbr_of_values; k++) {
        double  phiK = autoCorrelation(values, nbr_of_values, k);
        printf("phi_%d = %.7f, target = %.7f \t %.5f \n",k,phiK, target, d_abs(phiK - target));
        if(d_abs(phiK - target) <= threshold ) {
            s = k;
            break;
        }       
    }
    return s;
}

double d_abs(double d) {
    if(d < 0)
        d = -d;
    return d;
}
