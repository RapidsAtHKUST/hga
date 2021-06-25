#pragma once

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <cstdio>
#include <cstdlib>

#include "align.h"

using namespace std;

#define NGPUS 8

inline void check(cudaError_t result, char const *const func, const char *const file, int const line) {
    if (result) {
        fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line, static_cast<unsigned int>(result),
                cudaGetErrorString(result), func);
        exit(EXIT_FAILURE);
    }
}

#define CUCHK(val) check((val), #val, __FILE__, __LINE__)

double gpu_graph_align(Align aln, int start_id, int end_id);
