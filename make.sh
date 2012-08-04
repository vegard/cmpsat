set -e
set -u

g++ -Wall -O2 -o cmpsat cmpsat.cpp -lalglib
g++ -std=c++0x -Wall -O2 -o compare compare.cc -lalglib
