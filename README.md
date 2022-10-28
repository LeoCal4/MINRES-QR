# MINRES-QR 

A MATLAB implementation of the MINRES algorithm, using QR to decompose the T matrix and solve the resulting system accordingly to find the solution at each iteration.

The only two functions that need to be called directly are:
- `generate_problem_matrices`: generates a random direct graph and the relative problem matrices;
- `minres_qr`: runs the algorithm on the specified A and b.

Each of the main function on all the project files are througly commented, hence for more information on any of them simply consult the relative initial explanation. 

To replicate the experiments found in the report:
- run the script `generate_experiments_graphs.m`
- move the generated files in a new folder called `matrices`
- run the script `minres_full_experiments.m`, un-commenting the desired version to run
- collect the desired results
