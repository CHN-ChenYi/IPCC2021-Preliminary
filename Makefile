CC = mpicxx
CCFLAGS = -std=c++11 -O3 -fopenmp -march=znver2
SOURCE = SLIC.cpp SLIC.h
CASE = 2

main : default


default : $(SOURCE)
	$(CC) $(CCFLAGS) -fprofile-use=SLIC.profile -DMYMPI -o SLIC $<

gen : $(SOURCE)
	rm -rf SLIC.profile
	$(CC) $(CCFLAGS) -fprofile-generate=SLIC.profile -DMYMPI -o SLIC $<
	sbatch script.slurm
	squeue
	@echo -e "\e[1;31mGenerating profile, please execute 'make' after the job is done...\e[0m"


slct : $(SOURCE)
	$(CC) $(CCFLAGS) -fprofile-use=SLIC.profile -DSLCT -DMYMPI -o SLIC $<

slct-gen : $(SOURCE)
	rm -rf SLIC.profile
	$(CC) $(CCFLAGS) -fprofile-generate=SLIC.profile -DSLCT -DMYMPI -o SLIC $<
	sbatch script.slurm $(CASE)
	squeue


prof : $(SOURCE)
	$(CC) $(CCFLAGS) -fprofile-use=SLIC.profile -DSLCT -DMYMPI -DPROF -o SLIC $<

prof-gen : $(SOURCE)
	rm -rf SLIC.profile
	$(CC) $(CCFLAGS) -fprofile-generate=SLIC.profile -DSLCT -DMYMPI -DPROF -o SLIC $<
	sbatch script.slurm $(CASE)
	squeue


run :
	sbatch script.slurm $(CASE)
	squeue

clean :
	rm ./SLIC
