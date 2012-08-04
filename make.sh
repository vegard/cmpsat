set -e
set -u

g++ -Wall -O2 -Isrc -o cmpsat src/cmpsat.cpp -lalglib
g++ -std=c++0x -Wall -O2 -Isrc -o compare src/compare.cc -lalglib
