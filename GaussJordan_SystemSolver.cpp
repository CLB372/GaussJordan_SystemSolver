// Created by Chris Bryant, Jan. 2019. CLB372@cornell.edu
//
// This program solves an N x N linear system of equations using Gauss-Jordan Elimination to reduce the
// coefficient matrix to row canonical form (also known as Reduced Row Echelon Form, or RREF for short).
// Gauss-Jordan Elimination is more efficient than the recursive Cramer's Rule solution (posted as a separate
// repository in this GitHub), especially for N > 10.
//
// Text file input must be in matrix form. For example, for the following 3 x 3 system of equations
// containing 3 variables x, y, and z:
//
//			ax + by + cz = d
//			fx + gy + hz = i
//			jx + ky + mz = n
//
// The input must look like the following (without any white space to the left of each row or'/' characters:
//
//			a,b,c,d
//			f,g,h,i
//			j,k,m,n
//
// Within each row the numbers must be separated by commas as shown above. Both integral and decimal
// values are allowed.
//
// Gauss-Jordan Elimination consists of the following procedure:
//	(1) If the matrix is already in RREF, then stop.
//	(2) Find the 1st column from the left with a nonzero entry ('a'), and move that row to the top.
//	(3) Multiply that row by 1/a to create a leading 1.
//	(4) Subtract multiples of that row from all the other rows make each entry above and below the leading 1 zero.
//	(5) Repeat the above steps for each row.
//
//	Steps 2, 3, and 4 are the 3 'elementary row operations' that are described in linear algebra texts.
//	They do not change the identity/solutions of the system.
//
// After receiving the text file input and performing Gauss-Jordan Elimination on the matrix,
// the solutions are simply plucked from the last column and displayed to the user.
//
// The input can have any number of rows and columns, so long as the number of rows = the number of columns - 1.
// If the input is not N x (Nx1) numbers, the program will skip to an error message and terminate.

#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

// function prototypes
vector<vector<double>> getFileInput(string filename);
vector<vector<double>> gaussJordanElimination(vector<vector<double>> x);
bool rowCanonicalForm(vector<vector<double>> x);
vector<double> solveSystem_PreReducedMatrix(vector<vector<double>> x);
void printMatrix(vector<vector<double>> x);

int main()
{

	vector<vector<double>> inputMatrix;
	string filename;
	vector <double> solutions;

	// get the file name that contains the input
	cout << "\n\n  Enter the file name containing the N x (N+1) matrix representing an N x N system\n"
		<< "  of equations (please include the .txt extension):\n"
		<< "                                                       ";
	cin >> filename;
	inputMatrix = getFileInput(filename);
	cout << "\n\n  ";

	// display the matrix of numbers to the user for confirmation purposes
	for (int i = 0; i < inputMatrix.size(); i++)
	{
		for (int j = 0; j < inputMatrix[i].size(); j++)
		{
			cout << inputMatrix[i][j] << " ";
		}
		cout << "\n  ";
	}

	if (inputMatrix.size() == 0)
	{
		cout << "     ERROR: There were zero rows of numbers in your text file.\n\n";
	}
	else
	{
		if (inputMatrix[0].size() != inputMatrix.size() + 1)
		{
			cout << "     ERROR: The provided matrix of numbers is not an N x (N+1) matrix.\n\n";
		}
		else
		{
			// perform Gauss-Jordan Elimination on the inputMatrix
			inputMatrix = gaussJordanElimination(inputMatrix);

			// display the result to the user
			cout << "\n\n     RESULT: \n";
			solutions = solveSystem_PreReducedMatrix(inputMatrix);
			for (int i = 0; i < solutions.size(); i++)
			{
				cout << "              var" << i + 1 << " = " << solutions[i] << "\n\n";
			}
		}
	}
	cin.ignore(); cout << "     Created by Chris Bryant, Jan. 2019. CLB372@cornell.edu\n";  cin.ignore();

	return 0;
}


vector<vector<double>> getFileInput(string filename)
// This function reads the file 'filename' (a *.txt file) and returns a matrix of the file's input values.
// The input must be in the format described in the comments at the top of this program.
{
	vector<vector<double>> inputMatrix;
	ifstream myFile;
	string stringInput;
	myFile.open(filename);

	// for every row of text in the file myFile...
	while (getline(myFile, stringInput))
	{
		stringstream linestream(stringInput);
		string individualNumber;

		// increase the number of rows in inputMatrix by 1
		inputMatrix.resize(inputMatrix.size() + 1);

		// for every individual number captured in the current row of the text file...
		while (getline(linestream, individualNumber, ','))
		{
			// increase the number of columns in inputMatrix by 1
			inputMatrix[inputMatrix.size() - 1].resize(inputMatrix[inputMatrix.size() - 1].size() + 1);

			// convert individualNumber from a string to a double, and input
			// the double value into the currentinputMatrix element
			stringstream cast(individualNumber);
			cast >> inputMatrix[inputMatrix.size() - 1][inputMatrix[inputMatrix.size() - 1].size() - 1];
		}
	}

	return inputMatrix;
}


vector<double> solveSystem_PreReducedMatrix(vector<vector<double>> x)
// Takes an input matrix 'x' THAT MUST ALREADY BE IN REDUCED ROW ECHELON FORM (RREF),
// scans the right-hand column for the solutions to the system of equations, and
// returns a vector of those solutions
{
	vector<double> solutions;

	for (int i = 0; i < x.size(); i++)
	{
		solutions.resize(i + 1);
		solutions[i] = x[i][x[i].size() - 1];
	}

	return solutions;
}


vector<vector<double>> gaussJordanElimination(vector<vector<double>> x)
// Takes an input matrix 'x' and returns the RREF version of that matrix
// by performing Gauss-Jordan elimination on a local copy of'x'
{
	vector<vector<double>> y;
	double tempNum;
	bool stopLoop;

	// copy x to y
	for (int i = 0; i < x.size(); i++)
	{
		y.resize(i + 1);
		for (int j = 0; j < x[i].size(); j++)
		{
			y[i].resize(j + 1);
			y[i][j] = x[i][j];
		}
	}


	// now perform Gauss-Jordan Elimination on y

	// (1) if y is already in RREF, then stop and return y (handled by the condition in the next 'for' loop)
	// otherwise, proceed to steps (2) thru (4)

	// (2) going column by column, find a row with a nonzero entry in that column,
	//	   and move it to the top (do not do this for the last column, though)

	// a goes column by column to determine which entry/submatrix we're currently considering...
	for (int a = 0; a < y.size() && !rowCanonicalForm(y); a++)
	{
		// if the top row already has a leading nonzero entry, then no action 
		// is necessary for step (2)
		if (y[a][a] == 0)
		{
			// i goes row by row...
			stopLoop = false;
			for (int i = a + 1; i < y.size() && !stopLoop; i++)
			{
				if (y[i][a] != 0)
				{
					// move row i to the top of the submatrix currently under consideration
					// by swapping it with row a
					for (int j = 0; j < y[i].size(); j++)
					{
						tempNum = y[a][j];
						y[a][j] = y[i][j];
						y[i][j] = tempNum;
					}
					stopLoop = true;
				}
			}
		}

		// (3) multiply row 'a' by the reciprocal of the leading nonzero coefficient
		//     to create a leading 1
		tempNum = y[a][a];
		for (int i = 0; i < y[a].size(); i++)
		{
			y[a][i] = y[a][i] * (1 / tempNum);
		}
		// printMatrix(y);

		// (4) subtract multiples of row 'a' from the rows below it to make each entry
		//     below the leading 1 zero
		for (int i = 0; i < y.size(); i++)
		{
			if (i != a)
			{
				// multiply row i by the negative reciprocal of y[i][a] and add it to row a,
				// replacing row i with the result
				if (y[i][a] != 0)
				{
					tempNum = -1 / y[i][a];
					for (int j = 0; j < y[i].size(); j++)
					{
						y[i][j] = (y[i][j] * tempNum) + y[a][j];
					}
				}
			}
			// printMatrix(y);
		}
	}


	// -- now every element in the matrix should be zero except the last column and the diagonal
	// -- finally, we'll multiply each row by the reciprocal of its diagonal value to make
	//    the diagonal all 1s again
	for (int i = 0; i < y.size() && !rowCanonicalForm(y); i++)
	{
		tempNum = 1 / y[i][i];
		for (int j = 0; j < y[i].size(); j++)
		{
			y[i][j] *= tempNum;
		}
	}

	// return the reduced matrix
	return y;
}


bool rowCanonicalForm(vector<vector<double>> x)
// This function returns TRUE if the matrix 'x' is in Reduced Row Echelon Form (RREF),
// also known as Row Canonical Form.
//
// This function returns FALSE if it's not in RREF or if the matrix 'x'
// is not an N x (N+1) matrix.
{

	// test to make sure that the matrix 'x' is a valid N x (N+1) matrix
	for (int i = 0; i < x.size(); i++)
	{
		if (i != 0)
		{
			if (x.size() + 1 != x[i].size())
			{
				return false;
			}
		}

		if (x[i].size() != x.size() + 1)
		{
			return false;
		}
	}

	// -- now, assured of a valid 'x' input matrix, check x to see whether it's in RREF
	// -- check that it's all 1s along the diagonal and all 0s in every other entry
	//		(except for the right-most column)
	//    -- if it isn't, then return FALSE
	for (int i = 0; i < x.size(); i++)
	{
		for (int j = 0; j < x[i].size() - 1; j++)
		{
			if (i == j)
			{
				if (x[i][j] != 1) { return false; }
			}
			else
			{
				if (x[i][j] != 0) { return false; }
			}
		}
	}

	// if the function has not yet returned FALSE, then 'x' is indeed
	// in RREF, and the function should return TRUE
	return true;
}


void printMatrix(vector<vector<double>> x)
// prints a matrix for debugging purposes
{
	cout << "\n          ***BEGIN PRINT MATRIX FUNCTION:\n";

	for (int i = 0; i < x.size(); i++)
	{
		cout << "                                               ";
		for (int j = 0; j < x[i].size(); j++)
		{
			cout << x[i][j] << " ";
		}
		cout << "\n";
	}
	cout << "            ***END OF MATRIX FUNCTION\n";
}