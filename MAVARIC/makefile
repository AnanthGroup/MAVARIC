INC = -I ./Potentials/ -I ./source/
CXX = icc -std=c++11 -O2 -Wall

#Target definition
LINK_TARGET = mavaric

#List of objects

DIR = ./source/

vpath %.cpp source
vpath %.h source

OBJS = $(addprefix $(DIR), MonteCarlo.o\
  Sampling.o\
  Dynamics.o\
  Hamiltonian.o\
  Forces.o\
  ABM_MV_RPMD.o\
  RK4_MV_RPMD.o\
  main.o)

all: $(LINK_TARGET)


$(LINK_TARGET) : $(OBJS)
	$(CXX) $(INC) -o $@ $^


source/%.o :  %.cpp
	$(CXX) $(INC) -c -o $@ $<

clean:
	rm source/*.o
