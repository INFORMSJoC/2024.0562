# data
- - - 

## raw
A raw dataset downloaded from **RevLib** (https://www.revlib.org), an online resource for benchmarks within the domain of reversible and quantum circuit design.
> **RevLib**: *Wille, Robert, et al. "RevLib: An online resource for reversible functions and reversible circuits." 38th International Symposium on Multiple Valued Logic (ismvl 2008). IEEE, 2008.*

Each raw file provides: 
* The associated metadata of the Boolean reversible function
* The output part of the truth table, partially masked with `-` to represent the unspecified bits for incompletely specified functions.

## pckl
Processed data extracted from the raw dataset and saved via the Python pickle type, `[file_name].p`.
Each raw data is parsed and preprocessed in a ready-to-use dictionary, providing helpful mappings between the computational basis states (CBSs) and model indices (e.g., state nodes, commodity, etc.) utilized in the proposed formulation.

## common
A consolidated lookup table that contains necessary data descriptions for the entire benchmark function. Preprocessed and independently saved for convenience.
* `F-SM-INFO.p`: Basic attributes describing each reversible Boolean function
* `F-TLPQcoeff.p`: Preprocessed parameters to represent the constraints associated with problem-specific properties 
