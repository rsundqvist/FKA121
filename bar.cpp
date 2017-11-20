#include <stdio.h>
#include <stdlib.h>

double* init_an_array(int sz);
double* init_an_array_WRONG(int sz); // Using this function will give a segfault!

int main (int argc, char *argv[])
{
	for (int i = 0; i < argc; ++i)
	{
		printf("argv[%d] = %s\n", i, argv[i]);
	}
	
	int sz = argc > 1 ? atoi(argv[1]) : 10;
	//double *a = init_an_array(sz);
	double *a = init_an_array_WRONG(sz);
	for (int i = 0; i < sz; ++i)
	{
		printf("a[%d] = %.5f\n", i, a[i]);
	}	

	free (a); a = NULL;

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

double* init_an_array_WRONG(int sz)
{
	double ans[sz]; // Stack allocation
	for (int i = 0; i < sz; ++i)
	{
		ans[i] = i;
	}
	return ans;
}

