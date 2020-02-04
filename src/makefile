
CC = g++

ODIR = ../obj
OPT = -std=c++11 -fopenmp -Ofast

OBJ = model.o utils.o brent.o 
LIBS = -llapack
BOOST = -Iboost

all: ../DOLmod_MCMC ../CanADM ../USADM


%.o: %.cpp %.hpp
	$(CC) -c $(OPT) -o $@ $< $(BOOST) 

mcmc_abc.o: mcmc_abc.cpp
	$(CC) -c $(OPT) -o $@ $< $(BOOST)

../DOLmod_MCMC: mcmc_abc.o $(OBJ)
	$(CC) -g $(OPT) -o $@ $< $(OBJ) $(LIBS) $(BOOST)

../CanADM: solveCanODE.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)

../USADM: solveUSODE.cpp
	$(CC) -g $(OPT) -o $@ $< $(LIBS) $(BOOST)

clean:
	rm -f  *.o *~ 


