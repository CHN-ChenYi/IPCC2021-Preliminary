# IPCC2021-Preliminary

Only supports MPI when world_size is 2. If not, turn off MYMPI flag to use the single node version.

只支持 2 节点的 MPI，如果环境不是如此的话，请关闭 MYMPI 这一宏定义运行单节点版本。

## Prerequisites 环境准备

The following command will load all the necessary modules including:
* gcc/9.1.0
* intel/20.4.3
* mpi/intel/20.4.3

下一命令将加载上述三个必要的模块

``` bash
source ./setenv.sh
```

## Compile and Running 编译及运行

### Default version 默认版本（即比赛使用的版本）

Before compiling, please copy the test case into the root directory of the project and change m_spcount to the specified number.

在编译之前，请将测试样例拷贝至项目根目录下并将代码中的 m_spcount 改为该样例的指定值。

``` bash
make
make run
```

### SLCT version 可根据命令行参数切换输入样例的版本（仅供测试使用）

``` bash
make slct
make run CASE=1
make run CASE=2 # or just 'make run'
make run CASE=3
```
