#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stat.h"
#include "vec3.h"

int main()
{
    double v1[3] = {1, 0, 1};
    double v2[3] = {1, 2, 3};
    double x[3] = {1, 0, 0};
    double y[3] = {0, 1, 0};
    
    double u[3];
    double c[3];
    double d[3];
    double s[3];
    
    copy(v1, c);
    unit(c, u);
    diff(v1, v2, d);
    sum(v1, v2, s);
    return 0;
}
