#include "../include/align.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include "../include/align_gpu.cuh"
#include "../include/debug.h"

using namespace std;

void Align::input_graph(string gname) {
    ifstream graphf(gname);
    graphf >> num_v >> num_e;

    inv = (int *)malloc(sizeof(int) * num_e);
    inoff = (int *)malloc(sizeof(int) * (num_v + 1));
    outv = (int *)malloc(sizeof(int) * num_e);
    outoff = (int *)malloc(sizeof(int) * (num_v + 1));
    ref_graph = (char *)malloc(sizeof(char) * num_v);

    for (int i = 0; i < num_e; i++) {
        graphf >> inv[i];
    }
    for (int i = 0; i < num_v + 1; i++) {
        graphf >> inoff[i];
    }
    for (int i = 0; i < num_e; i++) {
        graphf >> outv[i];
    }
    for (int i = 0; i < num_v + 1; i++) {
        graphf >> outoff[i];
    }
    for (int i = 0; i < num_v; i++) {
        graphf >> ref_graph[i];
        if (ref_graph[i] == 'A') {
            ref_graph[i] = 1;
        } else if (ref_graph[i] == 'T') {
            ref_graph[i] = 2;
        } else if (ref_graph[i] == 'C') {
            ref_graph[i] = 4;
        } else if (ref_graph[i] == 'G') {
            ref_graph[i] = 8;
        } else {
            ref_graph[i] = 0;
        }
    }

    debug(num_v, num_e);
}

void Align::input_reads(string rname) {
    ifstream readf(rname);
    string tmp_str;

    num_reads = max_read_len = 0;

    while (readf >> tmp_str) {
        readf >> tmp_str;
        ++num_reads;
        max_read_len = max(max_read_len, (int)tmp_str.length());
    }

    debug(num_reads, max_read_len);

    readf.clear();
    readf.seekg(0);

    reads = (char *)malloc(sizeof(char) * num_reads * max_read_len);
    read_len = (int *)malloc(sizeof(int) * num_reads);

    for (int read_id = 0; read_id < num_reads; ++read_id) {
        char *cur_read = reads + read_id * max_read_len;
        readf >> tmp_str;
        readf >> cur_read;

        read_len[read_id] = strlen(cur_read);

#ifdef REVCOM
        reverse(cur_read, cur_read + strlen(cur_read));
        for (int i = 0; i < strlen(cur_read); i++) {
            if (cur_read[i] == 'A') {
                cur_read[i] = 'T';
            } else if (cur_read[i] == 'T') {
                cur_read[i] = 'A';
            } else if (cur_read[i] == 'C') {
                cur_read[i] = 'G';
            } else if (cur_read[i] == 'G') {
                cur_read[i] = 'C';
            }
        }
#endif

        for (int i = 0; i < strlen(cur_read); i++) {
            if (cur_read[i] == 'A') {
                cur_read[i] = 1;
            } else if (cur_read[i] == 'T') {
                cur_read[i] = 2;
            } else if (cur_read[i] == 'C') {
                cur_read[i] = 4;
            } else if (cur_read[i] == 'G') {
                cur_read[i] = 8;
            } else {
                cur_read[i] = 0;
            }
        }
    }
}

int Align::ind(int id) { return inoff[id + 1] - inoff[id]; }

int Align::outd(int id) { return outoff[id + 1] - outoff[id]; }

void Align::print_results() {
    for (int read_id = 0; read_id < num_reads; read_id++) {
        cout << read_id << " " << res_score[read_id] << " " << res_row[read_id] << " " << res_col[read_id] << endl;
    }
}

void Align::dump_graph() {
    cout << num_v << " " << num_e << endl;
    for (int i = 0; i < num_e; i++) {
        cout << inv[i] << endl;
    }
    for (int i = 0; i < num_v + 1; i++) {
        cout << inoff[i] << endl;
    }
    for (int i = 0; i < num_e; i++) {
        cout << outv[i] << endl;
    }
    for (int i = 0; i < num_v + 1; i++) {
        cout << outoff[i] << endl;
    }
    for (int i = 0; i < num_v; i++) {
        cout << ref_graph[i] << endl;
    }
}

void graph_align(Align aln) {
    aln.res_score = (int *)malloc(sizeof(int) * aln.num_reads);
    aln.res_row = (int *)malloc(sizeof(int) * aln.num_reads);
    aln.res_col = (int *)malloc(sizeof(int) * aln.num_reads);

    double wall_t = 1;

    if (aln.mode == 1) {  // GPU alignment
        wall_t = gpu_graph_align(aln, 0, aln.num_reads);
    }

    // aln.print_results();

    cout << "Time: " << wall_t << "s, GCUPS: " << (double)aln.num_v * aln.max_read_len * aln.num_reads / wall_t / 1000000000
         << endl;
}
