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

double velocityCorr(double (*vel)[3], double (*prevVel)[3], int N) {
    int i;
    double sum = 0;
    for(i = 0; i < N; i++) {
        sum += vel[i][0]*prevVel[i][0] + vel[i][1]*prevVel[i][1] + vel[i][2]*prevVel[i][2];
    }
    return sum/N;
}

double msdCorr(double (*pos)[3], double (*prevPos)[3], int N) {
    int i;
    double sum = 0, x1, y1, z1, x0, y0, z0;
    for(i = 0; i < N; i++) {
        x1 = pos[i][0];
        y1 = pos[i][1];
        z1 = pos[i][2];

        x0 = prevPos[i][0];
        y0 = prevPos[i][1];
        z0 = prevPos[i][2];


        sum += (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    }
    return sum/N;
}