set -e
set -u

g++ -Wall -O2 -o cmpsat cmpsat.cpp normaldistr.cpp ap.cpp
