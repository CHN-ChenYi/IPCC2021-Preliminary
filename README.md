# IPCC2021-Preliminary

Only supports MPI when world_size is 2. If not, turn off MYMPI flag to use the single node version.

## Compile

### Default version

```
make
```

### SLCT version
```
make slct
```

## Run

### Default version

```
make run
```

### SLCT version
```
make run CASE=1
make run CASE=2 # or just 'make run'
make run CASE=3
```
