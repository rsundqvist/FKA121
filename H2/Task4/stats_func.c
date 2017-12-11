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

// autocorrelation support function
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

// 
int findS(double * values, int nbr_of_values) {
    int k;
    int s = -1, extra = nbr_of_values, found = 0;     
    
    FILE *file;
	file = fopen("phi_k.dat","w");
	
	
    double target = 0.1353352832;
    for(k = 0; k < nbr_of_values; k++) {
        double  phiK = autoCorrelation(values, nbr_of_values, k);
        fprintf(file, "%d \t %.4f \n", k, phiK);
		extra--;
        if(phiK < target && !found) {
            s = k;
            found = 1;
            extra = k*20; // Generate extra data so the graph doesn't cut off.
        }
        if(extra <= 0) {
            break;
        } 
    }	
	fclose(file);
    
    return s;
}

double d_abs(double d) {
    return d < 0 ? -d : d;
}

// Get block average for a single B
double blockAverageS(double * f, int sz, int B) {
    int j, Mb = sz/B; // M_b = number of blocks of size B
    double F[Mb];
    for(j = 0; j < Mb; j++) {
        F[j] = blockAverage(f, sz, j+1, B); // +1 for i \in [1, B]
    }
    return B*sampleVariance(F,nba)/sampleVariance(f, sz);
}

// Block average for a single (B, j)-pair
double blockAverage(double * f, int nbr_of_values, int j, int B) {
    int i;
    double sum = 0;
    for(i = 0; i < B; i++)
        sum += f[i + (j-1)*B];
    return sum/B;
}

double sampleVariance(double * values, int nbr_of_values) {
    int i;
    double mean = getMean(values, nbr_of_values);
    double sum = 0;
    for(i = 0; i < nbr_of_values; i++)
        sum += (values[i] - mean)*(values[i] - mean);
    return sum/(nbr_of_values - 1);
}
