set -e
set -u

g++ -Wall -O2 -o cmpsat cmpsat.cpp normaldistr.cpp ap.cpp
g++ -std=c++0x -Wall -O2 -o compare compare.cc normaldistr.cpp ap.cpp
