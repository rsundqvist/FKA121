#include<stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "stats_func.h"

void readValues(const char *url, double * values, int N);

int main() {
    int N = 1000000; // 10^6
    double values[N];
    readValues("MC.txt", values, N);

    int s1 = findS(values, N, 0);
    printf("s = %d\n", s1);
    
    
    int B;
    double s2; 
    FILE *file;
	file = fopen("block_average.dat","w");
    for (B = 2; B < 30000; B++) {
        s2 = blockAverageS(values, N, B);
		fprintf(file, "%d \t %.4f \n", B, s2);
		if (B%3000==0) printf("B = %d\n", B);
    }
    fclose(file);
    printf("Done!\n");
}   



void readValues(const char *url, double *values, int N) {
    FILE * f = fopen(url, "r");
    if(f == NULL)
        printf("Failed to open file: \"%s\"", url);
 
    int i; 
    for(i = 0; i < N; ++i) {
        fscanf(f, "%lf",&values[i]);
    }
    close(f);
}
