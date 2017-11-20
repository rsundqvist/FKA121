#include<stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


double mathFunction1(double x);
void mcIntegrate(double *ans, double (*fuctionPtr)(double), double x, double x2, int N, gsl_rng * q);


int main() {
    const gsl_rng_type *T;
    gsl_rng *q;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    q = gsl_rng_alloc(T);
    gsl_rng_set(q,time(NULL));

    double ans[2];
    double (*Ptr)(double) = mathFunction1;

    int N = 10, i = 1;
    while(i <= 10) {
        mcIntegrate(ans, mathFunction1, 0, 1, N, q);
        printf("N^%d: \t mu = %.8f, sigma = %.8f \n", i, ans[0], ans[1]);
        N *= 10;
        i++;
    }

    // free memory
	gsl_rng_free(q);
}

double mathFunction1(double x) {
    return x - x*x;
}

void mcIntegrate(double *ans, double (*functionPtr)(double), double x1, double x2, int N, gsl_rng * q) {

    double rx, y;
    double ysum = 0, ysqSum = 0;
    // Sample N points
    int i;    
    for(i = 0; i < N; i++) {
        rx = x1 + (x2-x1)*gsl_rng_uniform(q);
        y = functionPtr(rx);
        ysum += y;
        ysqSum += y*y;
    }
    double mu = ysum/N;
    double sigmasq = 1/N * ysqSum - mu*mu;
    if (sigmasq < 0) sigmasq = -sigmasq;

    ans[0] = mu;
    ans[1] = sqrt(sigmasq/N);
}





