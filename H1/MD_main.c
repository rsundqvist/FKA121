/*
 MD_main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"

#define N 256

/* Main program */
int main()
{
    double  pos[N][3]; // Positions
    double  f[N][3]; // Forces
    int Nc = 10; // #primitive calls in each direction
    double a0 = 1; // Lattice parameter

    double L = 1; // Length of supercell

    /*
     Code for generating a uniform random number between 0 and 1. srand should only
     be called once.
    */
    //srand(time(NULL));
    double random_value = 0;
    //random_value = (double) rand() / (double) RAND_MAX;
    printf("random_value = %.5f", 0.0);

    /*
     Descriptions of the different functions in the files initfcc.c and
     alpotential.c are listed below.
    */

    /* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of 
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */
    init_fcc(pos, Nc, a0); // 
    
    /* 
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the 
     supercell and N is the number of atoms.
    */
    double energy;
    energy = get_energy_AL(pos, L, N);
    printf("energy = %.5f", energy);
    
    /* 
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell 
     and N is the number of atoms.
    */
    double virial;
    virial = get_virial_AL(pos, L, N);
    printf("virial = %.5f", virial);
    
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the 
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of 
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
    get_forces_AL(f,pos, L, N);  

    return 0; // Prevent err printout
}
