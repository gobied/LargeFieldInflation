g++ -c bkg.cpp -std=gnu++11 -lgsl -lgslcblas
ar rcs libbkg.a bkg.o
g++ main.cpp -std=gnu++11 -o main.o -lgsl -lgslcblas -L. -lbkg

