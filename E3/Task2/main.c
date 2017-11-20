#include<stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265359

double mathFunction1(double x);
double mathFunction2(double x);
double mathFunction3(double x);
double mathFunctionSinPi(double x);
double mathFunctionArcsinPi(double x);
void mcIntegrate(double *ans, double (*fuctionPtr)(double), double x, double x2, int N, gsl_rng * q);
void mcIntegrateImportanceSampling(double *ans, double (*pFun)(double), double (*invFun)(double), double x1, double x2, int N, gsl_rng * q);
double transformationMethod(double (*invFun)(double), double x1, double x2, gsl_rng * q);


gsl_rng * init_rng(); // gsl rng create
void free_rng(gsl_rng *);// gsl rng create+delete


int main() {
	gsl_rng * q = init_rng();

    double ans[2]; // Return value from mcIntegrate
    int N = 10, i = 1;
    // Bounds of inverse cummulative function
    double x1 = mathFunction3(0);
    double x2 = mathFunction3(1);
    while(i <= 5) {
        mcIntegrateImportanceSampling(ans, mathFunction2, mathFunction3, x1, x2, N, q);
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

double mathFunction2(double x) {
    return mathFunction1(x)/mathFunctionSinPi(x);
}

double mathFunctionSinPi(double x) {
    return 0.5 * PI * sin(PI*x);
}

double mathFunction3(double x) {
    return acos(1 - 2*x)/PI;
}

void mcIntegrate(double *ans, double (*fun)(double), double x1, double x2, int N, gsl_rng * q) {

    double rx, y;
    double ysum = 0, ysqSum = 0;
    // Sample N points
    int i;    
    for(i = 0; i < N; i++) {
        rx = x1 + (x2-x1)*gsl_rng_uniform(q);
        y = fun(rx);
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

void mcIntegrateImportanceSampling(double *ans, double (*gFun)(double), double (*invFun)(double), double x1, double x2, int N, gsl_rng * q) {
    double rx, y;
    double ysum = 0, ysqSum = 0, testsum = 0;
    // Sample N points
    int i;    
    for(i = 0; i < N; i++) {
        rx = transformationMethod(invFun, x1,x2,q);
        testsum += rx;
        y = gFun(rx);
        ysum += y;
        ysqSum += y*y;
    }
    printf("%.5f \n",testsum/N);
    double mu = ysum/N;
    double sigmasq = 1/N * ysqSum - mu*mu;
    if (sigmasq < 0) sigmasq = -sigmasq;

    ans[0] = mu;
    ans[1] = sqrt(sigmasq/N);
}


// Generate a random value from a continous probability distribution
double transformationMethod(double (*invFun)(double), double x1, double x2, gsl_rng * q) {
    
    double rx = x1 + (x2-x1)*gsl_rng_uniform(q);
    //printf("rx \t %.5f \t %.5f \n",rx, invFun(rx));
    return invFun(rx);
}


