#
# Makefile
#
all:
	g++ main.cc -lm -lfftw3 -lgsl -lgslcblas 
#-lboost_thread-mt
