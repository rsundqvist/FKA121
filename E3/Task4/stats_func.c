#include<stats_func.h>
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
