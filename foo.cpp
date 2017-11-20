#include <gsl/gsl_rng.h> //gsl must be installed!
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Stage 1 - Compilation
// Create foo.o from foo.c. -c flag tells gcc not to link
//gcc -c -o foo.o foo.c 
// Stage 2 - Link
//lm, lgsl, lgslcblas are libraries. -O3 (not a zero!) is Optimization3. Wall all warning
// "create executable run from foo.o"
//gcc -o run foo.o -O3 -Wall -lm -lgsl -lgslcblas 
// Stage 3 - Run executable
// ./run #run program

gsl_rng * init_rng();
void free_rng(gsl_rng *);

int main (int argc, char *argv[])
{
	int count = argc > 1 ? atoi(argv[1]) : 10;

	gsl_rng * q = init_rng();

	for (int i = 0; i < count; ++i)
	{
		double u = gsl_rng_uniform(q); // [0,1) - Note half-open intervall
		printf("%.5f\n",u);
	}

	printf("Goodbye!\n");
	return 0;
}

double* init_an_array(int sz)
{
	double *ans = malloc (sizeof(double) * sz); // Heap allocation
	for (int i = 0; i < sz; ++i)
	{
		ans[i] = i;
	}
	return ans;
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

void free_rng(gsl_rng * q)
{
	gsl_rng_free(q); // Free memory
}

