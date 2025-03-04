#ifndef PAPI_H
#define PAPI_H

#include <papi.h>

#include <stdexcept>

class PAPI {
public:
    static int setup(int operation_type);
    static long long int close(int event_set);
    static void start(int event_set);
};

#endif    // PAPI_H
