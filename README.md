# IPCC2021-Preliminary

​	目前用时
```plain
Computing time=1750 ms
```

​	目前程序的主要瓶颈在PerformSuperpixelSegmentation_VariableSandM函数，其他可以先不用管

​	我目前加了带线程锁的omp



Update:

​	目前用时700~800ms

```
[sca2301@ln121%bscc-a5 IPCC2021-Preliminary]$ make prof
g++ -std=c++11 -O3 -fopenmp -mavx2 -DPROF -o SLIC SLIC.cpp

[sca2301@ln121%bscc-a5 IPCC2021-Preliminary]$ make run
srun -p amd_256 -N 1 ./SLIC
Part1 time=360 ms
Part2 time=25 ms
Part3 time=23 ms
PerformSuperpixelSegmentation_VariableSandM time=463 ms
Computing time=730 ms
There are 0 points' labels are different from original file.
```

​	我把PerformSuperpixelSegmentation_VariableSandM函数分为Part1，Part2，Part3三部分。

​	Part2和Part3我加了OMP的Reduce快了很多，但是Part1加了OMP后（现在的是不用线程锁的，用线程锁的版本被ifdef了）效果不是很明显，从原先的3000+ms到360ms左右。Part计算部分有点奇怪，还需想一想怎么进一步优化。
