# HGA (Heterogeneous Graph Aligner)

This repository contains the source code of the paper "Accelerating Sequence-to-Graph Alignment on Heterogeneous Processors", by [Zonghao Feng](http://www.cse.ust.hk/~zfengah/) and Prof. [Qiong Luo](http://www.cse.ust.hk/~luo/).

## Build

Execute the following command to compile the project:

```
# Compile
make all

# Print debug information
make debug
```

## Run

Parameters of the program:
```
./bin/hga [-grmnodbt]
-g: graph file
-r: read file
-m: match score
-n: mismatch score
-o: gap penalty
-d: number of GPUs
-b: number of thread blocks
-t: number of threads per block
```

The suggested parameters for RTX 2080 Ti GPU are 68 thread blocks and 128 threads per block.

Scalability evaluation:
```
./bin/scale.sh 
```

## Ongoing work

- Extend to GPU clusters to support larger datasets
- Support more genome graph formats
- Intra-sequence parallelization
