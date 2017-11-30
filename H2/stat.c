#include <math.h>

double get_mean(double *vector, int sz) {
    double sum = 0;
    int j;
    for (j = 0; j < sz; ++j)
        sum+=vector[j];
    
    return sum/sz;
}

double get_variance(double *vector, double mean, int sz) {
    double sum = 0, d;
    int j;
    for (j = 0; j < sz; ++j) {
        d = vector[j] - mean;
        sum += d*d;
    }
    
    return sum/sz;
}

double get_mean_norm(double (*vector)[3], int sz) {
    double sum = 0, x, y, z;
    int j;
    for (j = 0; j < sz; ++j) {
        x = vector[j][0];
        y = vector[j][1];
        z = vector[j][2];
        sum += sqrt(x*x + y*y + z*z);
    }
    
    return sum/sz;
}

double get_variance_norm(double (*vector)[3], double mean, int sz) {
    double sum = 0, d, x, y, z;
    int j;
    for (j = 0; j < sz; ++j) {
        x = vector[j][0];
        y = vector[j][1];
        z = vector[j][2];
        d = sqrt(x*x + y*y + z*z) - mean;
        sum += d*d;
    }
    
    return sum/sz;
}
