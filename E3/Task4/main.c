#include<stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void readValues(const char *url, double * values, int N);

int main() {
    int N = 1000000; // 10^6
    double values[N];
    readValues("MC.txt", values, N);
    
    unsigned int i;
    for(i = 0; i < 15; i++)
        printf("%lf\n", values[i]);
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
