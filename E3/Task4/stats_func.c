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
    int s = -1, extra = nbr_of_values, found = 0;     
    
    FILE *file;
	file = fopen("phi_k.dat","w");
	
	
    double target = 0.1353352832;
    for(k = 0; k < nbr_of_values; k++) {
        double  phiK = autoCorrelation(values, nbr_of_values, k);
        //printf("phi_%d = %.7f, target = %.7f \t %.5f \n",k,phiK, target, phiK - target);
		fprintf(file, "%d \t %.4f \n", k, phiK);
		extra--;
        if(phiK < target && !found) {
            s = k;
            found = 1;
            printf("FOUND!\nFOUND!\nFOUND!\nFOUND!\nFOUND!\nFOUND!\nFOUND!\nFOUND!\nFOUND!\nFOUND!\n");
            extra = k*20;
        }
        if(extra <= 0) {
            break;
        } 
    }	
	fclose(file);
    
    return s;
}

double d_abs(double d) {
    if(d < 0)
        d = -d;
    return d;
}

double blockAverageS(double * values, int nbr_of_values, int B) {
    int k, j;
    int nBlockAverages = (int)round(nbr_of_values/B);
    double blockAverages[nBlockAverages];
    setZero(blockAverages, nBlockAverages);
    for(k = 0; k < nBlockAverages; k++) {
        j = k+1;
        blockAverages[k] = blockAverage(values, nbr_of_values, j, B);
    }
    return B *sampleVariance(blockAverages,nBlockAverages)/sampleVariance(values, nbr_of_values);
}

double blockAverage(double * values, int nbr_of_values, int j, int B) {
    int i;
    double sum = 0;
    for(i = 0; i < B; i++)
        sum += values[i + (j-1)*B];
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

void setZero(double * arr, int N) {
    int i;
    for (i = 0; i < N; i++){
        arr[i] = 0;
    }
}
