#include "vec3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double norm(double vec[3]) {
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

void copy(double vec[3], double copy[3]) {
    copy[0] = vec[0];
    copy[1] = vec[1];
    copy[2] = vec[2];
}

void diff(double vec1[3], double vec2[3], double diff[3]) {
    diff[0] = vec1[0] - vec2[0];
    diff[1] = vec1[1] - vec2[1];
    diff[2] = vec1[2] - vec2[2];
}

void sum(double vec1[3], double vec2[3], double sum[3]) {
    sum[0] = vec1[0] + vec2[0];
    sum[1] = vec1[1] + vec2[1];
    sum[2] = vec1[2] + vec2[2];
}

void unit(double vec[3], double unit[3]) {
    double d = norm(vec);
    
    unit[0] = vec[0]/d;
    unit[1] = vec[1]/d;
    unit[2] = vec[2]/d;
}

double dot(double vec1[3], double vec2[3]) {
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

double angle(double vec1[3], double vec2[3]) {
    double n = norm(vec1)*norm(vec2);
    return acos(dot(vec1, vec2)/n);
}

void print_vec(double vec[3]) {
    printf("[%.5f, %.5f, %.5f]", vec[0], vec[1], vec[2]);
}
