CC = mpicxx
CCFLAGS = -std=c++11 -O3 -fopenmp -mavx2
CASE = 1

main : SLIC.cpp SLIC.h
	$(CC) $(CCFLAGS) -o SLIC $<

prof : SLIC.cpp SLIC.h
	$(CC) $(CCFLAGS) -DSLCT -DPROF -o SLIC $<

slct : SLIC.cpp SLIC.h
	$(CC) $(CCFLAGS) -DSLCT -o SLIC $<

prof-lock : SLIC.cpp SLIC.h
	$(CC) $(CCFLAGS) -DPROF -DLOCK -o SLIC $<

gprof : SLIC.cpp SLIC.h
	$(CC) $(CCFLAGS) -pg -o SLIC $<

run :
	sbatch script.slurm --export=CASE=$(CASE)
	squeue

clean :
	rm ./SLIC
