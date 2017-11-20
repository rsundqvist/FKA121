#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

char* concat(const char *s1, const char *s2) {
	char *ans = malloc(strlen(s1)+strlen(s2)+1); // +1 for \0-terminator
	strcpy(ans, s1);
    strcat(ans, s2);
    return ans;
}

int periodic_bound(int i, int N) {
    if (i < 0) i += (-i/N+1) * N;
    return i % N;
}

double periodic_boundT(double i, double lim) {
    int tmpi = (int)(0.5+i*100000000);
    int tmplim = (int)(0.5+lim*100000000);
    if (tmpi < 0) tmpi += (-tmpi/tmplim+1) * tmplim;
    return (double)(tmpi % tmplim)/(100000000);
}
