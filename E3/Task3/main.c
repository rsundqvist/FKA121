#include<stdio.h>
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Exakt svar: 0.875

//Metropolis
void metropolisStep(double prev[3], double next[3], double Delta, double (*probFunc)(double[3],double[3]), gsl_rng * q);
double propFunction(double c1[3], double c2[3]);
double nextCoordinate(double c, double Delta, gsl_rng * q);
void generateMarkovChain(double (*chain)[3], double (*probFunc)(double[3],double[3]), double Delta, int N, gsl_rng * q);
void metropolisIntegrate3D(double *ans, double (*gfun)(double[3]), double(*samplePoints)[3], int N);
double mathFunction4(double r[3]);

void setZero(double (*arr)[3], int N);

// gsl stuff
gsl_rng * init_rng(); // gsl rng create



int count,total;
int main() {
	gsl_rng * q = init_rng();
	/*
	int j;
	for (j = 0; j < 1000; j++)
	    printf("r = %.5f\n", gsl_rng_uniform(q));

    return;*/
    double ans[2]; // Return value from integration
    int N = 10, i = 1;
    double Delta = 1.8;
    while(i <= 5) {
        count = 0, total = 0;
        double markovChain[N][3];
        //setZero(markovChain, N);
        generateMarkovChain(markovChain, propFunction, Delta, N, q);
        metropolisIntegrate3D(ans, mathFunction4, markovChain, N);
        printf("N=10^%d \t In = %.5f \t sigma = %.5f \n",i,ans[0],ans[1]);
        N *= 10;
        i++;
        //printf("%d    %d    %.3f \n", count, total, (double)count/total); // print accept rate. Target ~ 0.5
    }

    // free memory
	gsl_rng_free(q);
}

void setZero(double (*arr)[3], int N) {
    int i;
    for (i = 0; i < N; i++){
        arr[i][0] = 0;
        arr[i][1] = 0;
        arr[i][2] = 0;
    }
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
    
	//printf("r = %.5f, cInFunction = %.5f\n", r, c);
	//printf("return value = %.5f\n", c + Delta*(r - 0.5));
	
    return c + Delta*(r - 0.5);
}

void generateMarkovChain(double (*chain)[3], double (*probFunc)(double[3],double[3]), double Delta, int N, gsl_rng * q) {
    int i;
    //Randomize start position in [-10,10]
    for(i=0; i<3; i++)
        chain[0][i] = -2+4*gsl_rng_uniform(q);
        
    for(i = 0; i < 3; i++) {
        double * vec = chain[i];
        printf("\t(%.5f, %.5f, %.5f)\n", vec[0], vec[1], vec[2]);
    }
    for(i = 1; i< N; i++) {
        metropolisStep(chain[i-1], chain[i], Delta, probFunc, q);
        //printf("i = %d\n", i);
    }
    printf("-------------------\n");
    for(i = 0; i < 7; i++) {
        double * vec = chain[i];
        printf("\t(%.5f, %.5f, %.5f)\n", vec[0], vec[1], vec[2]);
    }
}

void metropolisStep(double prev[3], double next[3], double Delta, double (*probFunc)(double[3],double[3]), gsl_rng * q) {
    int i;
    double tmp[3];
    for (i = 0; i < 3; ++i) {
        //printf("prev[i] = %.5f\n", prev[i]);
        tmp[i] = nextCoordinate(prev[i], Delta, q);
    }
    double pr = probFunc(prev, tmp);
    double r = gsl_rng_uniform(q);
    
    if (pr > r) {// Accept new step
        for (i = 0; i < 3; ++i)
            next[i] = tmp[i];
    } else {
        for (i = 0; i < 3; ++i)
            next[i] = prev[i];
    }
    
    if (pr >= r)        
        count++;
    total++;
    //printf("%.3f (%.6f,%.6f,%.6f) \n\n", pr,next[0],next[1],next[2]);
}

void metropolisIntegrate3D(double *ans, double (*gFun)(double[3]), double(*samplePoints)[3], int N) {
    double f;
    double fsum = 0, fsqSum = 0;

    int i;    
    for(i = 0; i < N; i++) {
        f = gFun(samplePoints[i]);
        fsum += f;
        fsqSum += f*f;
    }
    double mu = fsum/N;
    double sigmasq = 1/N * fsqSum - mu*mu;
    if (sigmasq < 0) sigmasq = -sigmasq;

    ans[0] = mu;
    ans[1] = sqrt(sigmasq/N);
}

double propFunction(double c1[3], double c2[3]){
    double c1Sq = -c1[0]*c1[0] - c1[1]*c1[1] - c1[2]*c1[2];
    double c2Sq = -c2[0]*c2[0] - c2[1]*c2[1] - c2[2]*c2[2];
    return exp(c2Sq)/exp(c1Sq);
}

double mathFunction4(double r[3]) {
    double x = r[0], y = r[1], z = r[2];
    return x*x*(1 + y*y + y*y*z*z);  //*exp(-x*x-y*y-z*z);
}
