DIR=/home/andreas/Dokumente/Studium/tmp4/acp2015
CC=g++
#LIBDIR=$(DIR)/Libraries
#TOOLSDIR=$(DIR)/tools
#INCS=-I$(LIBDIR)
#INCS+=-I$(TOOLSDIR)
#CFLAGS=-std=c++11 -I. -O3 $(INCS)
CFLAGS=-std=c++11 -I. -O3 -Wall -Wextra

# the translation rule:
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

metropolis: main_mc.o metropolisState.o
	$(CC) -o metropolis main_mc.o metropolisState.o $(CFLAGS)

tests: test_MCbinning.o
	$(CC) -o test_binning test_MCbinning.o $(CFLAGS)