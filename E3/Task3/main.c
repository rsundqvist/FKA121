#include<stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265359

//Metropolis
void metropolisStep(double c[3], double ans[3], double Delta, double (*probFunc)(double[3],double[3]), gsl_rng * q);
double propFunction(double c1[3], double c2[3]);
double nextCoordinate(double c, double Delta, gsl_rng * q);
void generateMarkovChain(double (*chain)[3], double (*probFunc)(double[3],double[3]), double Delta, int nbr_of_points, gsl_rng * q);
void metropolisIntegrate3D(double *ans, double (*gfun)(double[3]), double(*samplePoints)[3], int nbr_of_points);
double mathFunction4(double r[3]);

// gsl stuff
gsl_rng * init_rng(); // gsl rng create


int main() {
	gsl_rng * q = init_rng();

    double ans[2]; // Return value from integration
    int N = 10, i = 1;
    double Delta = 10;
    while(i <= 5) {
        double markovChain[N][3];
        generateMarkovChain(markovChain, propFunction, Delta, N, q);
        metropolisIntegrate3D(ans, mathFunction4, markovChain, N);
        printf("N=10^%d \t In = %.5f \t sigma = %.5f \n",i,ans[0],ans[1]);
        N *= 10;
        i++;
    }

    // free memory
	gsl_rng_free(q);
}

double mathFunction4(double r[3]) {
    double x = r[0], y = r[1], z = r[2];
    return x*x + x*x*y*y + x*x*y*y*z*z;
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

double nextCoordinate(double c, double Delta, gsl_rng * q) {
    double r = gsl_rng_uniform(q); // [0, 1) - n.b. half-open!
    return c + Delta*(r - 0.5);
}

void generateMarkovChain(double (*chain)[3], double (*probFunc)(double[3],double[3]), double Delta, int nbr_of_points, gsl_rng * q) {
    int i;
    //Randomize start position in [-10,10]
    for(i=0; i<3; i++)
        chain[0][i] = gsl_rng_uniform(q);
    for(i = 1; i< nbr_of_points; i++)
        metropolisStep(chain[i],chain[i+1] , Delta, probFunc, q);
}

void metropolisStep(double c[3], double ans[3], double Delta, double (*probFunc)(double[3],double[3]), gsl_rng * q) {
    int i;
    double cNext[3];
    for (i = 0; i < 3; ++i) {
        cNext[i] = nextCoordinate(c[i], Delta, q);
    }
    double pr = probFunc(c, cNext);
    double r = gsl_rng_uniform(q);
    if (pr >= r) // Accept new step
        for (i = 0; i < 3; ++i)
            ans[i] = cNext[i];
    //printf("%.3f (%.6f,%.6f,%.6f) \n", pr,cNext[0],cNext[1],cNext[2]);
}

void metropolisIntegrate3D(double *ans, double (*gFun)(double[3]), double(*samplePoints)[3], int nbr_of_points) {
    double y;
    double ysum = 0, ysqSum = 0;

    int i;    
    for(i = 0; i < nbr_of_points; i++) {
        y = gFun(samplePoints[i]);
        ysum += y;
        ysqSum += y*y;
    }
    double mu = ysum/nbr_of_points;
    double sigmasq = 1/nbr_of_points * ysqSum - mu*mu;
    if (sigmasq < 0) sigmasq = -sigmasq;

    ans[0] = mu;
    ans[1] = sqrt(sigmasq/nbr_of_points);
}

double propFunction(double c1[3], double c2[3]){
    double c1Sq = -c1[0]*c1[0] - c1[1]*c1[1] - c1[2]*c1[2];
    double c2Sq = -c2[0]*c2[0] - c2[1]*c2[1] - c2[2]*c2[2];
    return exp(c2Sq)/exp(c1Sq);
}
