#
# Makefile
#
CC=g++
CFLAGS=
LDFLAGS= -lm -lfftw3 -lgsl -lgslcblas -fopenmp # -lfftw3_omp 

HDR= functs.hh params.hh cl_sim.hh
OBJ= main.o functs.o cl_sim.o

main.x : $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)

main.o : main.cc $(HDR)
	$(CC) $(CFLAGS) -c $<

functs.o : functs.cc params.hh cl_sim.hh
	$(CC) $(CFLAGS) -c $<

cl_sim.o : cl_sim.cc params.hh cl_sim.hh
	$(CC) $(CFLAGS) -c $<

clean :
	rm -f *.o *.x
