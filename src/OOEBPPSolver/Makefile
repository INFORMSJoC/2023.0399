###############################################################################################
######################### Compilação com suporte a múltiplos arquivos #########################
###############################################################################################

.PHONY: build clean

ARQ_MAIN = $(shell pwd | rev | cut -d '/' -f 1 | rev)

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
FLAGVERSION := $(shell gurobi_cl --version | cut -c 26,27,29 | head -n 1)
 
CPPFLAGS	=-g -O3 -Wall -Wextra -lm --std=c++17
INC      	= $(GUROBI_HOME)/include/
GUROBIFLAGS	=-I$(INC) -L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi$(FLAGVERSION)

build: $(OBJECTS)
	g++ *.o -o $(ARQ_MAIN) $(CPPFLAGS) $(GUROBIFLAGS)

%.o: %.cpp
	g++ -c $< -o $@ $(CPPFLAGS) -I$(INC)

softClean:
	rm $(ARQ_MAIN); rm -f *.o;

hardClean:
	rm $(ARQ_MAIN); rm -f *.o; rm -r out;
	