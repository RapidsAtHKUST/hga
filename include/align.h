#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Align {
   public:
    int match, mis, gap;
    int n_devices, grid_size, block_size, n_threads;
    int mode;

    int num_v, num_e;
    int *inv, *inoff, *outv, *outoff;
    char *ref_graph;

    int num_reads, max_read_len;
    char *reads;
    int *read_len;

    int *res_score, *res_row, *res_col;

    Align(int p_match = 2, int p_mis = 2, int p_gap = 3, int p_n_devices = 1, int p_grid_size = 68, int p_block_size = 128,
          int p_n_threads = 1, int p_mode = 1) {
        match = p_match;
        mis = p_mis;
        gap = p_gap;
        n_devices = p_n_devices;
        grid_size = p_grid_size;
        block_size = p_block_size;
        n_threads = p_n_threads;
        mode = p_mode;
    }

    int sub(char x, char y) { return (x == y) ? match : -mis; }

    void input_graph(string gname);

    void input_reads(string rname);

    int ind(int id);

    int outd(int id);

    void print_results();

    void dump_graph();
};

void graph_align(Align aln);
