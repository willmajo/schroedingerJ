***********************************
User manual for Schroedinger_solver
***********************************

1. Create an inputfile in the following format:
    2.0             #mass
    -2.0 2.0 1999   #xMin xMax nPoint
    1 5             #first and last eigenvalue to print
    linear          #interpolation type
    2               #nr. of interpolation points and xy declarations
    -2.0  0.0
     2.0  0.0

2. Save inputfile as schroedinger.inp in a new directory in the directory inputdata

3. Choose the directory and set the parameters for the illustration in schroedinger.py
    set options: * direc = directory of the inputfile
                 * bulge_factor = bulge amplitude of the wavefunctions
                 *lim_x = set the limits for the x-axis with additional summand
                 *lim_y = set the limits for the y-axis with additional summand

4. Execute ./schroedinger.py

5. View the illustations and outputfiles in the directory
