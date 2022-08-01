/* @file stat.h
**
** inline stat
** @@
******************************************************************************/

#ifndef STAT_H
#define STAT_H

#include <stdint.h>
#include <string.h>
#include <stdlib.h>

float ks_ksmall_float(size_t, float*, size_t);
int16_t ks_ksmall_int16_t(size_t, int16_t*, size_t);

static inline float meanf(const float *x, int n) {
    float sum = 0;
    int i = 0;
    for (i = 0; i < n; i++) {
        sum += x[i];
    }
    return sum / n;
}

static inline float meani16(const int16_t *x, int n) {
    float sum = 0;
    int i = 0;
    for (i = 0; i < n; i++) {
        sum += x[i];
    }
    return sum / n;
}


static inline float stdvf(const float *x, int n) {
    float sum = 0;
    int i = 0;
    float m = meanf(x, n); //can reuse already calculated mean
    for (i = 0; i < n; i++) {
        sum += (x[i] - m) * (x[i] - m);
    }
    return sqrtf(sum / n);
}

static inline float stdvi16(const int16_t *x, int n) {
    float sum = 0;
    int i = 0;
    float m = meani16(x, n); //can reuse already calculated mean
    for (i = 0; i < n; i++) {
        sum += (x[i] - m) * (x[i] - m);
    }
    return sqrtf(sum / n);
}

static inline float medianf(const float* x, int n) {

    float *copy = (float *)malloc(n * sizeof(float));
    memcpy(copy, x, n * sizeof(float));
    float m = ks_ksmall_float(n, copy, n / 2);
    free(copy);
    return m;

}

static inline int16_t mediani16(const int16_t* x, int n) {

    int16_t *copy = (int16_t *)malloc(n * sizeof(int16_t));
    memcpy(copy, x, n * sizeof(int16_t));
    int16_t m = ks_ksmall_int16_t(n, copy, n / 2);
    free(copy);
    return m;
}

#endif
