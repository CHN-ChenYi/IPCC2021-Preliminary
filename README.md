# IPCC2021-Preliminary

## Compile

### Default version

```
make gen # optimize by using runtime data on the current test case
make
```

### SLCT version
```
make slct-gen # optimize by using runtime data on the test case 2
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
