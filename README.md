# ShorAlgorithm
## Index
  - [Overview](#overview) 
  - [Setting](#setting)


## Overview
<!-- Write Overview about this project -->
- Shor algorithm for N6X5, N21X4, N35X16, N65X16
- Measurement error mitigation protocol available


## Features
**There are four class in the code: ShorN6X5, ShorN21X4, ShorN35X16, ShorN65X16**

### Setting
Need to setting initial account in ```BackendInit```
```
# backend setting
provider = IBMQ.load_account()
my_provider = IBMQ.get_provider("YOUR PROVIDER CODE")
backends = my_provider.backends()
```

### Basic Arguments
```
shots=1000, backend_name="ibmq_sasm_simulator, IsWtirl=False"
```

```SimulationQC()``` shows the probability distribution of the quantum phase estimation (order finding problem)
```
ShorN6X5(shots=1000, backend_name="ibm_nairobi", IsTwirl=False).SimulationQC()
```
NOTE: If ```IsTwirl=True```, then apply measurement error mitigation protocol

```RunShorAlgorithm()``` shows the result of the Shor's algorithm
```
ShorN6X5(shots=1000, backend_name="ibm_nairobi", IsTwirl=False).RunShorAlgorithm()
```
NOTE: This only available for ```IsTwirl=False```


### Setting Quantum Device
- ShorN6X5 is available for quantum device over 7-qubit
- ShorN21X4, ShorN35X16, ShorN65X16 are available for quantum device over 5-qubit
