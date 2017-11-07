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

int periodicBound (int i, int N) {
        if (i < 0) i += (-i/N+1) * N;
        return i % N;
}
