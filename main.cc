#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "include/align.h"
#include "include/debug.h"

using namespace std;

int main(int argc, char **argv) {
    int match = 2, mis = 2, gap = 3;
    int n_devices = 1, grid_size = 68, block_size = 128, n_threads = 1;
    int mode = 1;
    string gname, rname;

    int opt;
    while ((opt = getopt(argc, argv, "g:r:m:n:o:d:b:t:c:x:")) != -1) {
        switch (opt) {
            case 'g':
                gname = optarg;
                break;
            case 'r':
                rname = optarg;
                break;
            case 'm':
                match = atoi(optarg);
                break;
            case 'n':
                mis = atoi(optarg);
                break;
            case 'o':
                gap = atoi(optarg);
                break;
            case 'd':
                n_devices = atoi(optarg);
                break;
            case 'b':
                grid_size = atoi(optarg);
                break;
            case 't':
                block_size = atoi(optarg);
                break;
            case 'c':
                n_threads = atoi(optarg);
                break;
            case 'x':
                mode = atoi(optarg);
                break;
            default:
                abort();
        }
    }

    assert(!gname.empty());
    assert(!rname.empty());

    Align aln(match, mis, gap, n_devices, grid_size, block_size, n_threads, mode);

    debug(gname, rname, aln.match, aln.mis, aln.gap);
    debug(aln.n_devices, aln.grid_size, aln.block_size, aln.n_threads, aln.mode);

    aln.input_graph(gname);
    aln.input_reads(rname);

    graph_align(aln);

    return 0;
}
