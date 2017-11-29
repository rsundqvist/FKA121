#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// Compute Self diffusion coefficient using values from velocity autocorrelation data
double mathFunction1(double x, double * f, double * t);
void readValues(const char *url, double (*values)[2], int N);
void mcIntegrate(double *ans, double (*fun)(double, double *, double *), double x1, double x2, int N, gsl_rng * q, double * f, double * t);

gsl_rng * init_rng(); // gsl rng create
void free_rng(gsl_rng *);// gsl rng create+delete


int main() {
    int i;
    // Load data
    int nElements = 50000;
    double data[nElements][2];	
    double fValues[nElements];
    double tValues[nElements];

    readValues("tv_973.dat", data, nElements); 
    for(i = 0; i < nElements; i++) {
        fValues[i] = data[i][1];
        tValues[i] = data[i][0];
    }    

    // gsl random
    gsl_rng * q = init_rng();

    double ans[2]; // Return value from mcIntegrate
    int N = 10;
    i = 1;
    while(i <= 10) {
        mcIntegrate(ans, mathFunction1, tValues[0], tValues[nElements], N, q, fValues, tValues);
        printf("N^%d: \t mu = %.8f, sigma = %.8f \n", i, ans[0], ans[1]);
        N *= 10;
        i++;
    }

    // free memory
	gsl_rng_free(q);
    return 0;
}

double mathFunction1(double x, double * f, double * t) {
    int index = 0;    
    while(t[index] < x) {
        index++;
    }
    return f[index];
}

void mcIntegrate(double *ans, double (*fun)(double, double *, double *), double x1, double x2, int N, gsl_rng * q, double * f, double * t) {

    double rx, y;
    double ysum = 0, ysqSum = 0;
    // Sample N points
    int i;    
    for(i = 0; i < N; i++) {
        rx = x1 + (x2-x1)*gsl_rng_uniform(q);
        y = fun(rx,f,t);
        ysum += y;
        ysqSum += y*y;
    }
    double mu = ysum/N;
    double sigmasq = 1/N * ysqSum - mu*mu;
    if (sigmasq < 0) sigmasq = -sigmasq;

    ans[0] = mu;
    ans[1] = sqrt(sigmasq/N);
}

gsl_rng * init_rng()
{
	const gsl_rng_type *T;
	gsl_rng *q;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	q = gsl_rng_alloc(T);
	gsl_rng_set(q,time(NULL));
	return q;
}

void readValues(const char *url, double (*values)[2], int N) {
    FILE * f = fopen(url, "r");
    if(f == NULL)
        printf("Failed to open file: \"%s\"", url);
 
    int i; 
    for(i = 0; i < N; ++i) {
        fscanf(f, "%lf,%lf",&values[i][0], &values[i][1]);
    }
    //close(f);
}
