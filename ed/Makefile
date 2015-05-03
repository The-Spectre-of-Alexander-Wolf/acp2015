DIR=/home/andreas/Dokumente/Studium/tmp4/acp2015
CC=g++
LIBDIR=$(DIR)/Libraries
TOOLSDIR=$(DIR)/tools
INCS=-I$(LIBDIR)
INCS+=-I$(TOOLSDIR)
INCS+=-I/home/andreas/arpack++/include
CFLAGS=-std=c++11 -I. -O3 $(INCS)

# the translation rule:
%.o: %.cpp
	$(CC) $(CFLAGS) `python-config --cflags` -c $< -o $@

BoseHubbard: main_BoseHubbard.o BoseHubbardState.o
	$(CC) -o BoseHubbard main_BoseHubbard.o BoseHubbardState.o $(INCS) $(CFLAGS) `python-config --ldflags`
