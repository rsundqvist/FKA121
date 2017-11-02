/* 	Program for use on part two on discrete Fourier transforms in E1 
	Created by Martin Gren on 2014-10-21
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fft_func.h"
#include "../../utils/utils.h"
#define PI 3.141592653589
#define n 250 /*number of timesteps*/

int main ()
{
	/* timestep dt */
	double dt = 0.1;

	/* Declare and set Iteration parameter, data arrays and frequency */
	double data[n];
	double freq[n];
	double powspec_data[n];
	double f_values[] = {0, 1, 2};
	double phi_values[] = {0, PI/2};
	char info_string[100];

	/* declare file variables */
	FILE *file1;
	FILE *file2;	

	int round = 0;
	for (int i = 0; i < 3; ++i) {
		double f = f_values[i];
		for (int j = 0; j < 2; ++j) {
			double phi = phi_values[j];
			sprintf(info_string, "_%.2f_%.2f.dat", f, phi);


			/* calculate data points */
			for (int k = 0; k < n; k++)
				data[i] = cos(2*PI*f*dt*k + phi);


			/*Save function values in file*/
			file1 = fopen(concat("dat/function",info_string),"w");
			for (int k = 0; k < n; k++)	{
				fprintf (file1,"%e \t %e \n", k*dt, data[k]);
			}
			fclose(file1);

			/* make FFT (powerspectrum) */
			powerspectrum(data, powspec_data, n);
			powerspectrum_shift(powspec_data,n);
			fft_freq_shift(freq, dt, n);

			/*Save powerspectrum data in file */
			file2 = fopen(concat("dat/powerspectrum",info_string),"w");
			for (int k = 0; k < n; k++)	{
				fprintf (file2,"%e \t %e\n", freq[k], powspec_data[k]);
			}

			fclose(file2);

			printf("round =  %d\n", ++round);
		}
	}


	return 0;
}

