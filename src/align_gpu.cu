#include <omp.h>
#include <unistd.h>

#include <cassert>
#include <chrono>
#include <cstdio>

#include "../include/align_gpu.cuh"
#include "../include/debug.h"

using namespace std;

// __constant__ int cd_ref_graph[];

__global__ void align_kernel_inter_naive(int *d_dp_prev, int *d_dp, int *inv, int *inoff, char *ref_graph, char *reads,
                                         int *read_len, int *res_score, int *res_row, int *res_col, int num_v, int max_read_len,
                                         int batch_size, int match, int mis, int gap) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int skip = blockDim.x * gridDim.x;

    int *dp_prev = d_dp_prev + tid * num_v;
    int *dp = d_dp + tid * num_v;

    for (int read_id = tid; read_id < batch_size; read_id += skip) {
        int max_score = 0, max_row = 0, max_col = 0;
        char *cur_read = reads + read_id * max_read_len;

        memset(dp_prev, 0, num_v * sizeof(int));

        for (int i = 0; i < read_len[read_id]; i++) {
#pragma unroll 4
            for (int j = 0; j < num_v; j++) {
                int cur_score = 0;
                int sub_score = (ref_graph[j] == cur_read[i]) ? match : -mis;

                cur_score = max(cur_score, max(dp_prev[j] - gap, sub_score));

                for (int k = inoff[j]; k < inoff[j + 1]; k++) {
                    // #ifdef DEBUG
                    //                     printf("%d %d %d %d\n", i, j, k, inv[k]);
                    // #endif
                    cur_score = max(cur_score, max(dp_prev[inv[k]] + sub_score, dp[inv[k]] - gap));
                }
                dp[j] = cur_score;
#ifdef DEBUG
                printf("%d %d %d %d\n", i, j, sub_score, dp[j]);
#endif
                if (max_score <= cur_score) {
                    max_score = cur_score;
                    max_row = i;
                    max_col = j;
                }
            }
            int *tmp = dp_prev;
            dp_prev = dp;
            dp = tmp;
        }

        res_score[read_id] = max_score;
        res_row[read_id] = max_row;
        res_col[read_id] = max_col;

        // printf("%d %d %d %d %d\n", read_id, read_len[read_id], res_score[read_id], res_row[read_id], res_col[read_id]);
    }
}

// inter-seq parallel, coalesced
__global__ void align_kernel_inter(int *d_dp_prev, int *d_dp, const int *__restrict__ inv, const int *__restrict__ inoff,
                                   const char *__restrict__ ref_graph, const char *__restrict__ reads,
                                   const int *__restrict__ read_len, int *res_score, int num_v, int max_read_len,
                                   int batch_size, int match, int mis, int gap, const uint32_t *__restrict__ flags) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int skip = blockDim.x * gridDim.x;

    uint32_t flag = 0;

    // int *dp_prev = d_dp_prev + tid * num_v;
    // int *dp = d_dp + tid * num_v;
    int *dp_prev = d_dp_prev;
    int *dp = d_dp;

    for (int read_id = tid; read_id < batch_size; read_id += skip) {
        int max_score = 0;
        // int max_row = 0, max_col = 0;

        // memset(dp_prev, 0, num_v * sizeof(int));
        for (int j = 0; j < num_v; j++) {
            // dp_prev[j] = 0;
            dp_prev[j * skip + tid] = 0;
        }

        int cur_len = read_len[read_id];
        for (int i = 0; i < cur_len; i++) {
            char cur_read = reads[read_id + i * batch_size];  // coalesced
            // if (i % 8 == 0) {
            //     cur_read = reads[];
            // }

            for (int j = 0; j < num_v; j++) {
                int offset = j * skip;
                int sub_score = (ref_graph[j] == cur_read) ? match : -mis;
                int cur_score = max(0, max(dp_prev[offset + tid] - gap, sub_score));

                if (j % 32 == 0) {
                    flag = flags[j / 32];
                    // profile = graph[cur_read * num_v + j / 32];
                    // printf("%d %u\n", j, flag);
                }

                // sub_score = (profile & 1) ? match : -mis;

                if (flag & 1) {
                    cur_score = max(cur_score, max(dp_prev[offset - skip + tid] + sub_score, dp[offset - skip + tid] - gap));
                } else {
                    for (int k = inoff[j]; k < inoff[j + 1]; k++) {
                        cur_score =
                            max(cur_score, max(dp_prev[inv[k] * skip + tid] + sub_score, dp[inv[k] * skip + tid] - gap));
                    }
                }
                flag >>= 1;

                dp[offset + tid] = cur_score;

                // if (max_score <= cur_score) {
                //     max_score = cur_score;
                //     max_row = i;
                //     max_col = j;
                // }
                max_score = max(max_score, cur_score);
            }
            int *tmp = dp_prev;
            dp_prev = dp;
            dp = tmp;
        }
        res_score[read_id] = max_score;
        // res_row[read_id] = max_row;
        // res_col[read_id] = max_col;
    }
}

// inter-seq parallel with shared memory

// #define SMEM_SIZE 8192

#define BHEIGHT 8
#define BWIDTH 16
#define BSIZE 128

__global__ void align_kernel_inter_shared_naive(int *d_dp_prev, int *d_dp, int *inv, int *inoff, char *ref_graph, char *reads,
                                                int *read_len, int *res_score, uint32_t *flags, int num_v, int max_read_len,
                                                int batch_size, int match, int mis, int gap) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int skip = blockDim.x * gridDim.x;

    extern __shared__ int smem[];
    int *smem_ptr = smem;
    int *sdp = smem_ptr + threadIdx.x * BSIZE;

    int *dp_prev = d_dp_prev;
    int *dp = d_dp;

    // uint32_t flag = 0;

    for (int read_id = tid; read_id < batch_size; read_id += skip) {
        for (int j = 0; j < num_v; j++) {
            dp_prev[j * skip + tid] = 0;
        }

        // for (int j = BSIZE - BHEIGHT; j < BSIZE; j++) {
        //     dp[j] = 0;
        // }

        int max_score = 0;
        int cur_len = read_len[read_id];

        for (int i = 0; i < cur_len; i += BHEIGHT) {
            for (int j = 0; j < num_v; j++) {
                int offset = (j % BWIDTH) * BHEIGHT;
                int cur_graph = ref_graph[j];

                // k = 0
                {
                    int sub_score = (cur_graph == reads[read_id + i * batch_size]) ? match : -mis;

                    sdp[offset] = max(0, max(sub_score, dp_prev[j * skip + tid] - gap));

                    for (int v = inoff[j]; v < inoff[j + 1]; v++) {
                        sdp[offset] =
                            max(sdp[offset], max(dp_prev[inv[v] * skip + tid] + sub_score, sdp[(inv[v] % BWIDTH) * BHEIGHT] -
                            gap));
                    }

                    max_score = max(max_score, sdp[offset]);
                }

                for (int k = 1; (k < BHEIGHT) && (i + k < cur_len); k++) {
                    int sub_score = (cur_graph == reads[read_id + (i + k) * batch_size]) ? match : -mis;

                    sdp[offset + k] = max(0, max(sub_score, sdp[offset + k - 1] - gap));

                    for (int v = inoff[j]; v < inoff[j + 1]; v++) {
                        sdp[offset + k] = max(sdp[offset + k], max(sdp[(inv[v] % BWIDTH) * BHEIGHT - 1] + sub_score,
                                                                 sdp[(inv[v] % BWIDTH) * BHEIGHT] - gap));
                    }

                    max_score = max(max_score, sdp[offset + k]);
                }

                dp[j * skip + tid] = sdp[offset + BHEIGHT - 1];
            }

            int *tmp = dp_prev;
            dp_prev = dp;
            dp = tmp;
        }

        res_score[read_id] = max_score;
    }
}

__global__ void align_kernel_inter_shared(uint8_t *d_dp_prev, uint8_t *d_dp, uint32_t *inv, uint32_t *flags,
                                          uint32_t *ref_graph, uint32_t *reads, int *res_score, int num_v, int read_len,
                                          int batch_size, int match, int mis, int gap) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    int skip = blockDim.x * gridDim.x;

    extern __shared__ int smem[];
    // int *sdp = smem + threadIdx.x * BSIZE;

    uint8_t *dp_prev = d_dp_prev;
    uint8_t *dp = d_dp;

    uint32_t cur_reads = 0;
    uint32_t cur_graph = 0;
    uint32_t cur_ind = 0;

    for (int read_id = tid; read_id < batch_size; read_id += skip) {
        for (int j = 0; j < num_v; j++) {
            dp_prev[j * skip + tid] = 0;  // coalesced
        }
        // memset(dp_prev, 0, num_v * skip * sizeof(uint8_t));

        int max_score = 0;

        for (int i = 0; i < read_len; i += BHEIGHT) {
            cur_reads = reads[read_id + (i / 8) * batch_size];
            int index = 0;

            for (int j = 0; j < num_v; j++) {
                int offset = (j % BWIDTH) * BHEIGHT;
                int prev_offset = offset ? offset - BHEIGHT : (BWIDTH - 1) * BHEIGHT;
                uint32_t cur_read = cur_reads;

                if (j % 8 == 0) {
                    cur_graph = ref_graph[j / 8];
                } else {
                    cur_graph >>= 4;
                }

                if (j % 32 == 0) {
                    cur_ind = flags[j / 32];
                } else {
                    cur_ind >>= 1;
                }

                for (int k = 0; (k < BHEIGHT) && (i + k < read_len); k++) {
                    int sub_score = (cur_graph & cur_read & 0xF) ? match : -mis;
                    cur_read >>= 4;

                    // int cur_score = k ? sdp[offset + k - 1] : dp_prev[j * skip + tid];
                    int cur_score = k ? smem[threadIdx.x + blockDim.x * (offset + k - 1)] : dp_prev[j * skip + tid];
                    cur_score = max(cur_score - gap, max(0, sub_score));

                    if (cur_ind & 1) {
                        // int leftup = k ? sdp[prev_offset + k - 1] : dp_prev[(j - 1) * skip + tid];
                        // cur_score = max(cur_score, max(leftup + sub_score, sdp[prev_offset + k] - gap));

                        int leftup = k ? smem[threadIdx.x + blockDim.x * (prev_offset + k - 1)] : dp_prev[(j - 1) * skip + tid];
                        cur_score =
                            max(cur_score, max(leftup + sub_score, smem[threadIdx.x + blockDim.x * (prev_offset + k)] - gap));

                        // int leftup = k ? sdp[prev_offset + k - 1] : 0;
                        // cur_score = max(cur_score, sdp[prev_offset + k] - gap);
                    } else {
                        int neighbor = inv[index];
                        while (neighbor & 0xFF) {
                            int id = j - (neighbor & 0xFF);
                            // int leftup = k ? sdp[(id % BWIDTH) * BHEIGHT + k - 1] : dp_prev[id * skip + tid];
                            // cur_score = max(cur_score, max(leftup + sub_score, sdp[(id % BWIDTH) * BHEIGHT + k] - gap));
                            int leftup = k ? smem[threadIdx.x + blockDim.x * ((id % BWIDTH) * BHEIGHT + k - 1)]
                                           : dp_prev[id * skip + tid];
                            cur_score = max(
                                cur_score,
                                max(leftup + sub_score, smem[threadIdx.x + blockDim.x * ((id % BWIDTH) * BHEIGHT + k)] - gap));
                            neighbor >>= 8;
                        }
                    }

                    max_score = max(max_score, cur_score);
                    // sdp[offset + k] = cur_score;
                    smem[threadIdx.x + blockDim.x * (offset + k)] = cur_score;
                }

                if (!(cur_ind & 1)) {
                    index++;
                }

                if (i + BHEIGHT < read_len) {
                    // dp[j * skip + tid] = sdp[offset + BHEIGHT - 1];
                    dp[j * skip + tid] = smem[threadIdx.x + blockDim.x * (offset + BHEIGHT - 1)];
                }
            }

            uint8_t *tmp = dp_prev;
            dp_prev = dp;
            dp = tmp;
        }

        res_score[read_id] = max_score;
    }
}

// TODO: intra-seq parallel
__global__ void align_kernel_intra(int *d_dp_prev, int *d_dp, int *inv, int *inoff, char *ref_graph, char *reads, int *read_len,
                                   int *res_score, int *res_row, int *res_col, int num_v, int max_read_len, int batch_size,
                                   int match, int mis, int gap) {}

double gpu_graph_align(Align aln, int start_id, int end_id) {
    int dev_cnt;
    CUCHK(cudaGetDeviceCount(&dev_cnt));
    assert(dev_cnt > 0 && aln.n_devices <= dev_cnt);

    int tot_reads = end_id - start_id;

    // int dev_threads = min(batch_size, aln.grid_size * aln.block_size);
    int dev_threads = aln.grid_size * aln.block_size;
    int batch_l[NGPUS], batch_r[NGPUS], batch_size[NGPUS];

    omp_set_num_threads(NGPUS);

    cudaStream_t stream[NGPUS];

    uint8_t *d_dp_prev[NGPUS], *d_dp[NGPUS];
    int *d_res_score[NGPUS];
    uint32_t *d_flags[NGPUS];
    uint32_t *d_inv_c[NGPUS];
    uint32_t *d_ref_graph_c[NGPUS];
    uint32_t *d_reads_c[NGPUS];
    // int *d_inv[NGPUS], *d_inoff[NGPUS];
    // char *d_ref_graph[NGPUS], *d_reads[NGPUS];
    // int *d_read_len[NGPUS];
    // int *d_res_row[NGPUS], *d_res_col[NGPUS];

    debug(tot_reads, dev_threads);

    // prepare GCSR graph
    uint32_t *flags = (uint32_t *)calloc(aln.num_v / 32 + 1, sizeof(uint32_t));
    uint32_t bit = 1;
    int v_cnt = 0;
    for (int i = 0; i < aln.num_v; i++) {
        if (i % 32 == 0) {
            bit = 1;
            // flags[i / 32] = 0;
        }
        if (aln.ind(i) == 1 && aln.inv[aln.inoff[i]] == (i - 1)) {
            flags[i / 32] |= bit;
        } else {
            v_cnt++;
        }
        bit <<= 1;
    }

    uint32_t *inv_c = (uint32_t *)calloc(v_cnt, sizeof(uint32_t));
    int index = 0;
    for (int i = 0; i < aln.num_v; i++) {
        if (!(aln.ind(i) == 1 && aln.inv[aln.inoff[i]] == (i - 1))) {
            for (int j = aln.inoff[i], k = 0; j < aln.inoff[i + 1]; j++) {
                inv_c[index] |= ((i - aln.inv[j]) << (k * 8));
            }
            index++;
        }
    }

    uint32_t *ref_graph_c = (uint32_t *)calloc(aln.num_v / 8 + 1, sizeof(uint32_t));
    for (int i = 0; i < aln.num_v; i++) {
        ref_graph_c[i / 8] |= (aln.ref_graph[i] << ((i % 8) * 4));
    }

#pragma omp parallel for
    for (int dev_id = 0; dev_id < aln.n_devices; dev_id++) {
        CUCHK(cudaSetDevice(dev_id));
        CUCHK(cudaStreamCreateWithFlags(&stream[dev_id], cudaStreamNonBlocking));

        batch_l[dev_id] = dev_id * tot_reads / aln.n_devices;
        batch_r[dev_id] = (dev_id + 1) * tot_reads / aln.n_devices;
        batch_size[dev_id] = batch_r[dev_id] - batch_l[dev_id];

        debug(batch_l[dev_id], batch_r[dev_id], batch_size[dev_id]);

        char *rreads = (char *)malloc(sizeof(char) * batch_size[dev_id] * aln.max_read_len);
        char *readss = aln.reads + (start_id + batch_l[dev_id]) * aln.max_read_len;
        for (int i = 0; i < batch_size[dev_id]; i++) {
            for (int j = 0; j < aln.max_read_len; j++) {
                rreads[i + j * batch_size[dev_id]] = readss[i * aln.max_read_len + j];
            }
        }

        uint32_t *reads_c = (uint32_t *)calloc(batch_size[dev_id] * aln.max_read_len / 8 + 1, sizeof(uint32_t));
        for (int i = 0; i < batch_size[dev_id] * aln.max_read_len; i++) {
            reads_c[i / 8] |= (rreads[i] << ((i % 8) * 4));
        }

        // dp
        CUCHK(cudaMalloc(&d_dp_prev[dev_id], sizeof(uint8_t) * aln.num_v * dev_threads));
        CUCHK(cudaMalloc(&d_dp[dev_id], sizeof(uint8_t) * aln.num_v * dev_threads));

        // graph
        // CUCHK(cudaMalloc(&d_inv[dev_id], sizeof(int) * aln.num_e));
        // CUCHK(cudaMalloc(&d_inoff[dev_id], sizeof(int) * (aln.num_v + 1)));
        // CUCHK(cudaMalloc(&d_ref_graph[dev_id], sizeof(char) * aln.num_v));
        CUCHK(cudaMalloc(&d_inv_c[dev_id], sizeof(uint32_t) * v_cnt));
        CUCHK(cudaMalloc(&d_ref_graph_c[dev_id], sizeof(uint32_t) * (aln.num_v / 8 + 1)));

        // read
        // CUCHK(cudaMalloc(&d_reads[dev_id], sizeof(char) * batch_size[dev_id] * aln.max_read_len));
        // CUCHK(cudaMalloc(&d_read_len[dev_id], sizeof(int) * batch_size[dev_id]));
        CUCHK(cudaMalloc(&d_reads_c[dev_id], sizeof(uint32_t) * (batch_size[dev_id] * aln.max_read_len / 8 + 1)));

        // result
        CUCHK(cudaMalloc(&d_res_score[dev_id], sizeof(int) * batch_size[dev_id]));
        // CUCHK(cudaMalloc(&d_res_row[dev_id], sizeof(int) * batch_size[dev_id]));
        // CUCHK(cudaMalloc(&d_res_col[dev_id], sizeof(int) * batch_size[dev_id]));

        CUCHK(cudaMalloc(&d_flags[dev_id], sizeof(uint32_t) * (aln.num_v / 32 + 1)));

        // CUCHK(cudaMemset(d_dp_prev, 0, sizeof(int) * aln.num_v * dev_threads));

        // CUCHK(cudaMemcpyAsync(d_inv[dev_id], aln.inv, sizeof(int) * aln.num_e, cudaMemcpyHostToDevice, stream[dev_id]));
        // CUCHK(cudaMemcpyAsync(d_inoff[dev_id], aln.inoff, sizeof(int) * (aln.num_v + 1), cudaMemcpyHostToDevice, stream[dev_id]));
        // CUCHK(cudaMemcpyAsync(d_ref_graph[dev_id], aln.ref_graph, sizeof(char) * aln.num_v, cudaMemcpyHostToDevice, stream[dev_id]));
        // CUCHK(cudaMemcpyAsync(d_reads[dev_id], aln.reads + (start_id + batch_l[dev_id]) * aln.max_read_len, sizeof(char) * batch_size[dev_id] * aln.max_read_len, cudaMemcpyHostToDevice, stream[dev_id]));
        // CUCHK(cudaMemcpyAsync(d_reads[dev_id], rreads, sizeof(char) * batch_size[dev_id] * aln.max_read_len, cudaMemcpyHostToDevice, stream[dev_id]));
        // CUCHK(cudaMemcpyAsync(d_read_len[dev_id], aln.read_len + (start_id + batch_l[dev_id]), sizeof(int) * batch_size[dev_id], cudaMemcpyHostToDevice, stream[dev_id]));

        CUCHK(cudaMemcpyAsync(d_reads_c[dev_id], reads_c, sizeof(uint32_t) * (batch_size[dev_id] * aln.max_read_len / 8 + 1),
                              cudaMemcpyHostToDevice, stream[dev_id]));
        CUCHK(cudaMemcpyAsync(d_flags[dev_id], flags, sizeof(uint32_t) * (aln.num_v / 32 + 1), cudaMemcpyHostToDevice,
                              stream[dev_id]));
        CUCHK(cudaMemcpyAsync(d_inv_c[dev_id], inv_c, sizeof(uint32_t) * v_cnt, cudaMemcpyHostToDevice, stream[dev_id]));
        CUCHK(cudaMemcpyAsync(d_ref_graph_c[dev_id], ref_graph_c, sizeof(uint32_t) * (aln.num_v / 8 + 1),
                              cudaMemcpyHostToDevice, stream[dev_id]));

        free(rreads);
        free(reads_c);
    }

    free(flags);
    free(inv_c);
    free(ref_graph_c);

    auto t_start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
    for (int dev_id = 0; dev_id < aln.n_devices; dev_id++) {
        CUCHK(cudaSetDevice(dev_id));

        // align_kernel_inter<<<aln.grid_size, aln.block_size>>>(
        //     d_dp_prev[dev_id], d_dp[dev_id], d_inv[dev_id], d_inoff[dev_id], d_ref_graph[dev_id], d_reads[dev_id],
        //     d_read_len[dev_id], d_res_score[dev_id], d_res_row[dev_id], d_res_col[dev_id], aln.num_v, aln.max_read_len,
        //     batch_size[dev_id], aln.match, aln.mis, aln.gap, d_flags[dev_id]);

        // align_kernel_inter<<<aln.grid_size, aln.block_size>>>(d_dp_prev[dev_id], d_dp[dev_id], d_inv[dev_id],
        // d_inoff[dev_id],
        //                                                       d_ref_graph[dev_id], d_reads[dev_id], d_read_len[dev_id],
        //                                                       d_res_score[dev_id], aln.num_v, aln.max_read_len,
        //                                                       batch_size[dev_id], aln.match, aln.mis, aln.gap,
        //                                                       d_flags[dev_id]);

        cudaFuncSetAttribute(align_kernel_inter_shared, cudaFuncAttributeMaxDynamicSharedMemorySize, 65536);
        // sizeof(int) * BSIZE * aln.block_size

        // align_kernel_inter_shared_naive<<<aln.grid_size, aln.block_size, 65536>>>(
        //     d_dp_prev[dev_id], d_dp[dev_id], d_inv[dev_id], d_inoff[dev_id], d_ref_graph[dev_id], d_reads[dev_id],
        //     d_read_len[dev_id], d_res_score[dev_id], d_flags[dev_id], aln.num_v, aln.max_read_len, batch_size[dev_id],
        //     aln.match, aln.mis, aln.gap);

        align_kernel_inter_shared<<<aln.grid_size, aln.block_size, 65536>>>(
            d_dp_prev[dev_id], d_dp[dev_id], d_inv_c[dev_id], d_flags[dev_id], d_ref_graph_c[dev_id], d_reads_c[dev_id],
            d_res_score[dev_id], aln.num_v, aln.max_read_len, batch_size[dev_id], aln.match, aln.mis, aln.gap);
    }

#pragma omp parallel for
    for (int dev_id = 0; dev_id < aln.n_devices; dev_id++) {
        CUCHK(cudaSetDevice(dev_id));

        CUCHK(cudaDeviceSynchronize());
    }

    auto t_end = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
    for (int dev_id = 0; dev_id < aln.n_devices; dev_id++) {
        CUCHK(cudaSetDevice(dev_id));

        CUCHK(cudaMemcpyAsync(aln.res_score + start_id + batch_l[dev_id], d_res_score[dev_id], sizeof(int) * batch_size[dev_id],
                              cudaMemcpyDeviceToHost, stream[dev_id]));
        // CUCHK(cudaMemcpyAsync(aln.res_row + start_id + batch_l[dev_id], d_res_row[dev_id], sizeof(int) * batch_size[dev_id],
        //                       cudaMemcpyDeviceToHost, stream[dev_id]));
        // CUCHK(cudaMemcpyAsync(aln.res_col + start_id + batch_l[dev_id], d_res_col[dev_id], sizeof(int) * batch_size[dev_id],
        //                       cudaMemcpyDeviceToHost, stream[dev_id]));
    }

#pragma omp parallel for
    for (int dev_id = 0; dev_id < aln.n_devices; dev_id++) {
        CUCHK(cudaSetDevice(dev_id));

        CUCHK(cudaStreamDestroy(stream[dev_id]));

        CUCHK(cudaFree(d_dp_prev[dev_id]));
        CUCHK(cudaFree(d_dp[dev_id]));
        // CUCHK(cudaFree(d_inv[dev_id]));
        // CUCHK(cudaFree(d_inoff[dev_id]));
        // CUCHK(cudaFree(d_ref_graph[dev_id]));
        // CUCHK(cudaFree(d_reads[dev_id]));
        // CUCHK(cudaFree(d_read_len[dev_id]));
        CUCHK(cudaFree(d_res_score[dev_id]));
        // CUCHK(cudaFree(d_res_row[dev_id]));
        // CUCHK(cudaFree(d_res_col[dev_id]));
        CUCHK(cudaFree(d_flags[dev_id]));
        CUCHK(cudaFree(d_inv_c[dev_id]));
        CUCHK(cudaFree(d_ref_graph_c[dev_id]));
        CUCHK(cudaFree(d_reads_c[dev_id]));
    }

    std::chrono::duration<double> diff = t_end - t_start;
    return diff.count();
}
