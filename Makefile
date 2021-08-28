CC = mpicxx
CCFLAGS = -std=c++11 -O3 -fopenmp -mavx2
SOURCE = SLIC.cpp SLIC.h
CASE = 1

main : $(SOURCE)
	$(CC) $(CCFLAGS) -o SLIC $<

prof : $(SOURCE)
	$(CC) $(CCFLAGS) -DSLCT -DPROF -o SLIC $<

mpi : $(SOURCE)
	$(CC) $(CCFLAGS) -DSLCT -DMYMPI -o SLIC $<

mpi-prof : $(SOURCE)
	$(CC) $(CCFLAGS) -DSLCT -DMYMPI -DPROF -o SLIC $<

slct : $(SOURCE)
	$(CC) $(CCFLAGS) -DSLCT -o SLIC $<

gprof : $(SOURCE)
	$(CC) $(CCFLAGS) -pg -o SLIC $<

run :
	sbatch script.slurm --export=CASE=$(CASE)
	squeue

clean :
	rm ./SLIC
