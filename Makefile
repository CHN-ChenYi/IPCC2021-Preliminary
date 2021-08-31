CC = mpicxx
CCFLAGS = -std=c++11 -O3 -fopenmp -march=znver2
SOURCE = SLIC.cpp SLIC.h
CASE = 2

main : default

default : $(SOURCE)
	$(CC) $(CCFLAGS) -DMYMPI -o SLIC $<

slct : $(SOURCE)
	$(CC) $(CCFLAGS) -DSLCT -DMYMPI -o SLIC $<

prof : $(SOURCE)
	$(CC) $(CCFLAGS) -DSLCT -DMYMPI -DPROF -o SLIC $<

run :
	sbatch script.slurm $(CASE)
	squeue

clean :
	rm ./SLIC
