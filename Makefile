CC=g++
NVCC=nvcc
CFLAGS=-O3# -DREVCOM
CUFLAGS=-arch=compute_75 -code=sm_75 -Xcompiler -fopenmp

PROG=hga
BIN=bin
INC=include
OBJS=$(BIN)/main.o $(BIN)/align.o $(BIN)/align_gpu.o

.PHONY:all debug clean

all: $(PROG)

$(PROG): $(OBJS)
	$(NVCC) $(CFLAGS) $(CUFLAGS) $(OBJS) -o $(BIN)/$@ 

$(BIN)/%.o: %.cc
	$(NVCC) $(CFLAGS) -c $< -o $@

$(BIN)/%.o: src/%.cc $(INC)/%.h
	$(NVCC) $(CFLAGS) $(CUFLAGS) -c $< -o $@

$(BIN)/%.o: src/%.cu $(INC)/%.cuh
	$(NVCC) $(CFLAGS) $(CUFLAGS) -c $< -o $@

debug: CFLAGS += -g -Xcompiler -Wall -DDEBUG
debug: all

clean:
	rm -f $(BIN)/*
