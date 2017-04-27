cc = g++
OPENMP = -fopenmp
CFLAGS = -03
LIBS = 


TARGETS = serial OPENMP

all: $(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o

