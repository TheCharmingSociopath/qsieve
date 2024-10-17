The files 'hashing_parameters_classical.py' and 'hashing_parameters_quantum.py' compute the number of hash tables t and filter parameter alpha for the dimensions of interest.
The file 'compute_data_heuristics.py' computes the number of physical qubits and time complexities for all sieves using the hash parameters previously computed.
The file 'make_plots.py' makes the plots based on the outputs of 'compute_data_heuristics.py'
The file 'auxiliary_functions.py' contens several auxiliary functions used to compute the hash parameters and number of physical qubits and time complexities
The remaining files are adapted from Litinski's repository in https://github.com/litinski/magicstates and concern magic state distillation protocols
All our raw data is in the folder 'data'
