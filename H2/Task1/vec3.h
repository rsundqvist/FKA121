
#ifndef _vec3_h
#define _vec3_h

double norm(double vec[3]);
double angle(double vec1[3], double vec2[3]);

void copy(double vec[3], double copy[3]);
void unit(double vec[3], double unit[3]);
void diff(double vec1[3], double vec2[3], double diff[3]);
void sum(double vec1[3], double vec2[3], double sum[3]);

double dot(double vec2[3], double vec1[3]);

void print_vec(double vec[3]);
#endif
