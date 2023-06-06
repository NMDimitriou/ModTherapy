# Makefile for Drug Resistance simulation model

IDIR:=include
ODIR:=obj
_DEPS:=declarations.hpp
DEPS:=$(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=ppp.o read_dat.o 
OBJ:=$(patsubst %,$(ODIR)/%,$(_OBJ))

MAIN:=ppp

all: $(MAIN)

CC=g++
RM=rm -rf

CXXFLAGS=-I$(IDIR) -O2 -static -std=c++17

$(ODIR)/%.o: %.cpp $(DEPS) $(DEPSP)
	$(CC) -c -o $@ $< $(CXXFLAGS) 

$(MAIN):$(OBJ)
	$(CC) -o $@ $^ $(CXXFLAGS) -Wall -Wextra -Werror

clean:
	$(RM) $(ODIR)/*.o *~ core*
clear:
	$(RM) *.txt ../DSR*
