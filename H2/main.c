#include <stdio.h>
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stat.h"
#include "vec3.h"

#define N 100
#define TIME_MAX 20000
#define PI 3.14159265359

int main()
{
    double v1[3] = {1, 0, 1};
    double v2[3] = {1, 2, 3};
    double v3[3] = {1, 0, 0};
    
    double u[3];
    double c[3];
    double d[3];
    double s[3];
    
    copy(v1, c);
    unit(c, u);
    diff(v1, v2, d);
    sum(v1, v2, s);
    
    printf("v1 = %s\n", to_string(v1));
    printf("v2 = %s\n", to_string(v2));
    printf("v3 = %s\n", to_string(v3));
    printf("norm(v2) = %.5f\n", norm(v2));
    printf("unit = %s\n", to_string(u));
    printf("copy = %s\n", to_string(c));
    printf("dot(v1, v2) = %.5f\n", dot(v1, v2));
    printf("diff(v1, v2) = %s\n", to_string(d));
    printf("sum(v1, v2) = %s\n", to_string(s));
    return 0;
}
