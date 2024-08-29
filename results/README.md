# results
- - - 
This directory contains a full set of results that are reported in the paper.
Each subdirectory corresponds to a different type of result files, while the file name shares a common naming format including the information of each experiment:
* Function name: [*func_name*]
* Default gate count: [*gate_count*]
* Extra gate count: [*extra_gate_count*]
* Configuration ID: [*config_ID*]

The table below presents the correspondence of configuration ID between the repository and the paper.

| Repository | E1 | E2 | E3 | E4 | E5 | E6 | E1T  | E2T  | E3T  | E4T  | E5T  | E6T  |
|------------|----|----|----|----|----|----|------|------|------|------|------|------|
| Paper      | C1 | C2 | C3 | C4 | C5 | C6 | C1BP | C2BP | C3BP | C4BP | C5BP | C6BP |

With the configuration ID a naming format of the result files in each subdirectory:

* sol: `[func_name]-D[gate_count]+[extra_gate_count]-E[config_ID].sol`
* log: `[func_name]-D[gate_count]+[extra_gate_count]-E[config_ID].txt`
* circuit: `[func_name]-D[gate_count]+[extra_gate_count]-E[config_ID].jpg`

## sol
This subdirectory contains the result files that save the final solutions retrieved from the Gurobi solver.
The files can be read by the ```pickle``` package, where the value of each single decision variable is saved as a nested dictionary 
in the form of ```{var_name : {var_index : var_value}}```. 

Each subdictionary is keyed with the alphabet (possibly greek) assigned to each variable type, 
For the indexed variables, they are keyed again with the tuple of associated indices so that the value can contain the final value of the corresponding variable.

> ```{'Z': 14.0, 'x': {('S', '000', 1, 0): 1.0, ('S', '000', 2, 0): 0.0, ... ```

## log
This subdirectory contains `.txt` files for all the logs that printed out from Gurobi solver 
at the time when the final experiment was conducted.

## circuit
This subdirectory contains `.jpg` files that represent the solutions as the corresponding circuit diagrams.