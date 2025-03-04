#include "PAPI.h"

int PAPI::setup(int operation_type) {
    if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
        throw std::runtime_error("PAPI Library initialization error!");
    }

    int event_set = PAPI_NULL;
    if(PAPI_create_eventset(&event_set) != PAPI_OK ||
       PAPI_add_event(event_set, operation_type) != PAPI_OK) {
        throw std::runtime_error("PAPI EventSet creation error!");
    }

    return event_set;
}

long long int PAPI::close(int event_set) {
    long long int fp_ops;
    if(PAPI_stop(event_set, &fp_ops) != PAPI_OK) {
        throw std::runtime_error("PAPI Stop error!");
    }

    if(PAPI_cleanup_eventset(event_set) != PAPI_OK ||
       PAPI_destroy_eventset(&event_set) != PAPI_OK) {
        throw std::runtime_error("PAPI Cleanup/Destroy error!");
    }

    return fp_ops;
}

void PAPI::start(int event_set) {
    if(PAPI_start(event_set) != PAPI_OK) {
        throw std::runtime_error("PAPI Start error!");
    }
}
