WITH_OPENMP ?= 0

ifeq ($(WITH_OPENMP),1)
  CXX = /opt/homebrew/opt/llvm/bin/clang++
  CXXFLAGS = -O3 -march=native -std=c++17 -fopenmp
else
  CXX = /usr/bin/clang++
  CXXFLAGS = -O3 -march=native -std=c++17
endif

all: chain_solver

chain_solver: chain_solver.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

interface_certify: interface_certify.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f chain_solver interface_certify

test: chain_solver
	./chain_solver 127
	./chain_solver 8191

verify: chain_solver
	./chain_solver --verify

.PHONY: all clean test verify
