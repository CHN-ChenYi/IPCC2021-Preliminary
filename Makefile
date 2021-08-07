CC = g++
CCFLAGS = -std=c++11 -O3 -fopenmp

main : SLIC.cpp SLIC.h
	$(CC) $(CCFLAGS) -o SLIC $<

run :
	srun -p amd_256 -N 1 ./SLIC

clean :
	rm ./SLIC
