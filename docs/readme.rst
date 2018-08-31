***********************************
User manual for Schroedinger_solver
***********************************

1. Create an inputfile in the following format:
    * 2.0               #mass
    * -2.0 2.0 1999     #xMin xMax nPoint
    * 1 5               #first and last eigenvalue to print
    * linear            #interpolation type
    * 2                 #nr. of interpolation points and xy declarations
    * -2.0  0.0
    *  2.0  0.0

2. Move to the schroedinger_project directory in your Unix shell.

3. Save inputfile as schroedinger.inp in a new directory in the directory inputdata.

4. Execute ./solver.py to solve the Schroedinger equation. Choose the directory with your 
   inputfile with the -d argument. The resualts will be saved in that directory as .dat files.

5. Execute ./visualize.py to illustrate the results. Choose the directory with your results with the  -d argument. The illustration will be saved in the directory as solution_schroedinger.pdf. For a better illustration, set the following options as floats:
    * -bf = --bulge_factor: bulge amplitude of the wavefunctions
    * -xl = --lim_x: set the limits for the x-axis with additional summand
    * -yl = --lim_y: set the limits for the y-axis with additional summand
