This program solves an N x N linear system of equations using Gauss-Jordan Elimination to reduce the
coefficient matrix to row canonical form (also known as Reduced Row Echelon Form, or RREF for short).
Gauss-Jordan Elimination is more efficient than the recursive Cramer's Rule solution (posted as a separate
repository in this GitHub), especially for N > 10.

The GaussJordan_SystemSolver.exe is the executable file for this program. A text file containing the input augmented matrix must be in the same folder as this file.

The GaussJordan_SystemSolver.cpp is the relevant C++ code.

Text file input must be in augmented matrix form. For example, for the following 3 x 3 system of equations
containing 3 variables x, y, and z:

			ax + by + cz = d
			fx + gy + hz = i
			jx + ky + mz = n

The input must look like the following (without any white space to the left of each row):

			a,b,c,d
			f,g,h,i
			j,k,m,n

Within each row the numbers must be separated by commas as shown above. Both integral and decimal
values are allowed.

Gauss-Jordan Elimination consists of the following procedure:
(1) If the matrix is already in RREF, then stop.
(2) Find the 1st column from the left with a nonzero entry ('a'), and move that row to the top.
(3) Multiply that row by 1/a to create a leading 1.
(4) Subtract multiples of that row from all the other rows make each entry above and below the leading 1 zero.
(5) Repeat the above steps for each row.

Steps 2, 3, and 4 are the 3 'elementary row operations' that are described in linear algebra texts.
They do not change the identity/solutions of the system.

After receiving the text file input and performing Gauss-Jordan Elimination on the matrix,
the solutions are simply plucked from the last column and displayed to the user.

The input can have any number of rows and columns, so long as the number of rows = the number of columns - 1.
If the input is not N x (Nx1) numbers, the program will skip to an error message and terminate.
