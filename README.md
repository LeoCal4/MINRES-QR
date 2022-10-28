# MINRES-QR 

A MATLAB implementation of the MINRES algorithm, using QR to decompose the T matrix and solving the resulting system accordingly to find the solution at each iteration.

The only two functions that need to be called directly are:
- `generate_problem_matrices`: generates a random direct graph and the relative problem matrices (truncated);
- `minres_qr`: runs the algorithm on the specified A and b.

Each of the main function on all the project files are througly commented, hence for more information on any of them simply consult the relative initial explanation. 


The `generate_experiments_graphs.m` and `minres_experiments.m` scripts have been employed to conduct the experiments shown in the report, but the Python functions used to calculate time/iterations averages and to generate plots from the experiments results are not included, as they are poorly written for the most part. Either way, the two scripts are maintained, as they may be useful for further experimentation or results verification.
