# 2024.0562


This repository provides a snapshot of source codes that are utilized to produce the computational results presented in **Jung and Choi (2024)**.

> (_In press_) Jung, Jihye and In-Chan Choi. “A New Multicommodity Network Flow Model and Branch-and- Cut for Optimal Quantum Boolean Circuit Synthesis.” INFORMS Journal on Computing.



## Dataset

### Raw Dataset
A raw dataset is saved in directory ```../data/raw/``` as `.txt` format.
The entire set of files are also available for download on **RevLib** (https://www.revlib.org), an online resource for benchmarks within the domain of reversible and quantum circuit design.
> **RevLib**: *Wille, Robert, et al. "RevLib: An online resource for reversible functions and reversible circuits." 38th International Symposium on Multiple Valued Logic (ismvl 2008). IEEE, 2008.*

### Processed Dataset
We extracted necessary information from the raw text files and saved the processed data as pickle files ```../data/processed/{function_name}.p```.

This dataset enables the user to retrieve a ready-to-use dictionary containing helpful mappings between the computational basis states (CBSs) and model indices (e.g., state nodes, commodity, etc.) used in the formulation.

### Common Info
For user convenience, we saved consolidated lookup tables under ```../data/common``` to store necessary function descriptions for the entire benchmark.
* `func_attrs.p`: Basic attributes describing each reversible Boolean function.
* `func_coeffs.p`: Preprocessed parameters representing the constraints associated with problem-specific properties utilized in the proposed formulation.


## Replications
Below is a brief description of how to replicate the computational results presented in the paper.

For experiment identification, we denote each attributes as follows:
* Function name: `{func_name}`
* Default gate count: `{def_d}`
* Extra gate count: `{extra_d}`
* Configuration ID: `{config_id}`

For detailed information about `{config_ID}`, refer to the table below presenting the correspondence of configuration ID between the repository and the paper.

| Repository | E1 | E2 | E3 | E4 | E5 | E6 | E1T  | E2T  | E3T  | E4T  | E5T  | E6T  |
|------------|----|----|----|----|----|----|------|------|------|------|------|------|
| Paper      | C1 | C2 | C3 | C4 | C5 | C6 | C1BP | C2BP | C3BP | C4BP | C5BP | C6BP |


### Requirements
* Python 3.8+
* Gurobi Optimizer 9.1.2 (install `gurobipy` for Python API and ensure you have a usable license)

### Model & Modeler Usage
We already generated the optimization models in `.mps` format in `../models`.
Below is the naming protocal used to identify each model.
```
{model_type}-{function_name}-D{default_gates}+{extra_gates}.mps
```
Particularly, for `{model_type}`, we offer two options: `MW` and `MS`. `MW` refers to a surrogate model while `MS` implies a strong model with entire _xi-cuts_ handled as original constraints. 

You can also generate a new mps model by using `../src/modeler.py` by the command below:
```
python src/modeler.py --func="{function_name}" --model_type="{model_type}" --extra_d={int}
```
> NOTE: This command will save a corresponding mps file in `../models` with a **TEST** label added at the end, such as `{file_name}-TEST.mps`.


### Solver Usage
This command below runs our Gurobi-based solver by retrieving the assigned `.mps` model from `..\models`.
```
python src/solver.py --func="{function_name}" --extra_d={int} --env_id="{config_ID}"
```
A single experiment will be conducted with the assigned experiment settings, saving the optimal solution and its log file into corresponding directory, `../results/sol` and `../results/log`.
> NOTE: Currently, this command will save a corresponding mps file in `../results` with a **TEST** label added at the end, such as `{file_name}-TEST.{extension}`.



## Results

This directory contains the entire set of result files that are reported in our paper.
Each subdirectory corresponds to a different type of result files, while the file name shares a common naming format including the information of each experiment:

**Solution File:** `../results/sol`
```
[func_name]-D[gate_count]+[extra_gate_count]-E[config_ID].sol`
```
This `pickle` file contains the final variable values retrieved from the Gurobi solver.
The values are saved in a nested dictionary format `{var_name : {var_index : var_value}}`. 

Each subdictionary is first keyed with the alphabet assigned to each variable type.
Then, the dictionary is keyed again with the tuple of associated indices representing the corresponding decision variable.

```
{'Z': 14.0, 'x': {('S', '000', 1, 0): 1.0, ('S', '000', 2, 0): 0.0, ...
```

**Log File:** `../results/log`
```
[func_name]-D[gate_count]+[extra_gate_count]-E[config_ID].txt
```

**Circuit Diagram:** `../results/circuit`
```
[func_name]-D[gate_count]+[extra_gate_count]-E[config_ID].jpg
```
This subdirectory contains `.jpg` files that represent the solutions as the corresponding circuit diagrams.