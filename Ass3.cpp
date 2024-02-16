/*
File:Ass3.cpp
Author: Jason Waid | 400424500
McMaster University | Numerical Linear Algebra and Numerical Optimization | SFWRTECH 4MA3
Professor: SESHASAI SRINIVASAN
Description: Implementation of the Broyden's Method
*/

#include <iostream>
#include <ios>
#include <iomanip>  
#include <limits>
#include <string>
#include <algorithm>
#include <iterator>
#include <utility>

using namespace std;

void performGaussianElimination(double matrix[][3], int numUnknowns);
void performBackwardsSubstitution(double matrix[][3], int matrixSize);
double calcTwoNorm(double matrix[2], int vectorSize);

/*
Function: displayBroydensMethodTableHeaders
Purpose: prints table headers for Rayleigh Quotient Iteration With Shifts
Param: NONE
Return: NONE
Author: Jason Waid
Date Modified: 03/26/2022
*/
void displayBroydensMethodTableHeaders() {
	cout << setw(2) << "k" << setw(25) << "x[k]^T" << setw(30) <<"B"<< setw(30)<< "f(x[k])"<< endl;
	cout << setfill('-') << setw(120) << "-" << endl;
	cout << setfill(' ');
}

/*
Function: performGaussianElimination
Purpose: Performs Gaussian Elimination on the matrix
Param: double matrix[][11] which is the matrix in question, int numRows which are the number of rows, int numColumns which are the number of columns
Return: NONE
Author: Jason Waid
Date Modified: 03/13/2022
*/
void performGaussianElimination(double matrix[][3], int numUnknowns) {
	double ratio = 0;

	for (int k = 0; k < numUnknowns - 1; k++) {
		if (matrix[k][k] == 0) {
			//stop
			exit(0);
		}
		for (int j = k + 1; j < numUnknowns; j++) {
			ratio = matrix[j][k] / matrix[k][k];

			for (int i = k; i < numUnknowns + 1; i++) {
				matrix[j][i] = matrix[j][i] - (ratio * matrix[k][i]);
			}

		}
	}
}

/*
Function: performBackwardsSubstitution
Purpose: Performs backwards substituion on the matrix
Param: double matrix[][] which is the matrix in question, int numUnknowns is the number of unknowns
Return: NONE
Author: Jason Waid
Date Modified: 03/16/2022
*/
void performBackwardsSubstitution(double matrix[][3], int matrixSize) {

	double temp = 0;
	double unknownVals[2] = { 0 };

	for (int i = 0; i < matrixSize; i++) {
		unknownVals[i] = matrix[i][matrixSize];
	}

	for (int i = matrixSize - 1; i >= 0; i--)
	{

		temp = 0;
		for (int j = i + 1; j < matrixSize; j++) {

			temp += matrix[i][j] * unknownVals[j];

		}

		if (matrix[i][i] == 0) {
			unknownVals[i] = 0;
		}
		else {
			unknownVals[i] = (matrix[i][matrixSize] - temp) / matrix[i][i];
		}

	}

	for (int i = 0; i < matrixSize; i++) {
		matrix[i][matrixSize] = unknownVals[i];
	}
}

/*
Function: calcTwoNorm
Purpose: calculated two norm of given vector
Param: double matrix: the given matrix, b values are given vector, matrixSize: size of matrix, also the column number of the b vals
Return: NONE
Author: Jason Waid
Date Modified: 03/14/2022
*/
double calcTwoNorm(double vector[2], int vectorSize) {

	double norm = 0;

	//sumation of all vals squared
	for (int i = 0; i < vectorSize; i++) {

		norm = norm + (fabs(vector[i]) * fabs(vector[i]));

	}

	//power of 1/2
	norm = sqrt(norm);


	return norm;
}


int main()
{
	string inputBuffer;
	int matrixSize = 0;
	int menu = 0;
	//double matrix[10][10 + 1] = {};
	double matrix[10][10 + 1] = { {1, 2, -2, 0}, {1, 4, -4, 0} };
	//							{ {x1.1,x1.2,const,0},{x2.1,x2.2,const,0}}
	double jacobian[10][10] = { {1, 0}, {0, 1} };

	//following the example using the same assumptions
	double Yk_transpose[2] = { 2, 3 };
	double Sk_transpose[2] = { 1, 1 };
	double Xk_transpose[2] = { 1, 2 };

	double new_Xk_transpose[2] = { 0,0 };

	double temp_vector_transpose[2] = {0};
	double temp_matrix[2][3] = { {0, 0, 0}, {0, 0, 0} };
	double temp_Fx_vector_result[2] = { 0,0 };
	double temp_NEW_Fx_vector_result[2] = { 0,0 };
	double temp_Val = 0;

	//top and botton of B[k+1] calculation
	double top_Val = 0;
	double bottom_Val = 0;

	matrixSize = 3;

	int k = 0;

	displayBroydensMethodTableHeaders();
	while (true) {

		//print details
		cout << setw(2) << k;

		cout << setprecision(4) << fixed << setw(10) <<"[ " << Xk_transpose[0] << ", " << Xk_transpose[1] << " ]";
		cout << "\t| ";
		cout << setw(10) << "[";

		for (int i = 0; i < 2; i++) {

			cout << "[ ";

			for (int j = 0; j < 2; j++) {

				cout << jacobian[i][j];
				if (j == 0) {
					cout << ", ";
				}	
			}

			cout << " ]";

			if (i == 0) {
				cout << ", ";
			}
		}

		cout << "]\t";

		
		//cout << "Yk^t: [ " << Yk_transpose[0] << "\t" << Yk_transpose[1] << " ]\n";

		//Step 1
		//solve BkSk = -f(Xk) for Sk using Gauss Elimination

		//part 1 = solve -f(Xk)
		//x + 2y - 2 = 0
		//x2^2 + 4y2^2 - 4 = 0

		//i think xk is x and y // yes this is correct

		
		cout << "[ ";

		for (int i = 0; i < 2; i++) {
			//reset val
			temp_Fx_vector_result[i] = 0;

			if (i == 0) {

				for (int j = 0; j < 3; j++) {

					if (j == 2) {

						temp_Fx_vector_result[i] += matrix[i][j];
					}
					else {
					
					temp_Fx_vector_result[i] += matrix[i][j] * Xk_transpose[j];

					}
				}

				cout << temp_Fx_vector_result[i] << ", ";

				//temp_Fx_vector_result[i] *= -1;
			}
			else {

				for (int j = 0; j < 3; j++) {

					if (j == 2) {

						temp_Fx_vector_result[i] += matrix[i][j];
					}
					else {
						temp_Fx_vector_result[i] += matrix[i][j] * (Xk_transpose[j] * Xk_transpose[j]);

					}
				}

				cout << temp_Fx_vector_result[i];
				//temp_Fx_vector_result[i] *= -1;
			}
		}

		cout << " ]\n";

		//pass expected val was [-3, -13]
		//cout << temp_Fx_vector_result[0] << " " << temp_Fx_vector_result[1] << "\n";

		//create temp matrix by combining Bk and temp_Fx_vector_result to pass into Gauss elimination to find Sk
		for (int row = 0; row < 2; row++) {
			for (int col = 0; col < 3; col++) {

				if (col == 2) {
					temp_matrix[row][col] = -temp_Fx_vector_result[row];
					//cout << "row: " << row << " col: " << col << "\n";
					//cout << temp_matrix[row][col];
				}
				else {
					temp_matrix[row][col] = jacobian[row][col];
				}
			}
		}


		//PASS, gave expected matrix
		//for (int row = 0; row < 2; row++) {
		//	cout << "[ ";
		//	for (int col = 0; col < 3; col++) {

		//		
		//		cout << temp_matrix[row][col] << " ";
		//		
		//	}
		//	cout << " ]\n";
		//}

		//Use Gauss elimination to find Sk
		performGaussianElimination(temp_matrix, 2);
		performBackwardsSubstitution(temp_matrix, 2);

		//resulting Sk should be [-3, -13]
		Sk_transpose[0] = temp_matrix[0][2];
		Sk_transpose[1] = temp_matrix[1][2];
		//PASS
		//cout << "Sk = [ " << Sk_transpose[0] << " " << Sk_transpose[1] << " ]\n";

		//Step 2
		//Xk+1 = Xk + Sk
		for (int row = 0; row < 2; row++) {

			new_Xk_transpose[row] = Xk_transpose[row] + Sk_transpose[row];

		}

		//Step 3
		//Yk = f(Xk+1) - f(Xk)

		//calc f(Xk+1)
		for (int i = 0; i < 2; i++) {
			//reset val
			temp_NEW_Fx_vector_result[i] = 0;

			if (i == 0) {

				for (int j = 0; j < 3; j++) {

					if (j == 2) {
						temp_NEW_Fx_vector_result[i] += matrix[i][j];
					}
					else{

						temp_NEW_Fx_vector_result[i] += matrix[i][j] * new_Xk_transpose[j];

					}
				}

			}
			else {

				for (int j = 0; j < 3; j++) {

					if (j == 2) {
						temp_NEW_Fx_vector_result[i] += matrix[i][j];
					}else{
						temp_NEW_Fx_vector_result[i] += matrix[i][j] * (new_Xk_transpose[j] * new_Xk_transpose[j]);

					}
				}

			}
		}

		//calc f(Xk)
		//temp_Fx_vector_result[]


		for (int row = 0; row < 2; row++) {

			Yk_transpose[row] = temp_NEW_Fx_vector_result[row] - temp_Fx_vector_result[row];

		}

		//Step 4
		//(Yk - Bk*Sk), we'll do Bk*Sk first

		//Bk * Sk
		for (int row = 0; row < 2; row++) {
			temp_vector_transpose[row] = 0;
			for (int col = 0; col < 2; col++) {

				temp_vector_transpose[row] += (jacobian[row][col] * Sk_transpose[row]);

			}
		}

		//PASS
		/*cout << "\nRESULT OF Bk * Sk:\n";
		cout << temp_vector_transpose[0] << "\t" << temp_vector_transpose[1] << "\n";*/

		//above results in vector

		//next, Yk - above
		//continue using temp val for calculation
		for (int row = 0; row < 2; row++) {

			temp_vector_transpose[row] = (Yk_transpose[row] - temp_vector_transpose[row]);

		}
		//above changes vector val

		//PASS
	/*	cout << "\nRESULT OF Yk - (Bk * Sk):\n";
		cout << temp_vector_transpose[0] << "\t" << temp_vector_transpose[1] << "\n";*/

		//top_Val is so far temp_vector_transpose[] - actual transpose of Sk (1x2)

		bottom_Val = 0;
		//calc bottom_Val AKA: 1/bottom_Val
		for (int i = 0; i < 2; i++) {
			
			bottom_Val += (Sk_transpose[i] * Sk_transpose[i]);

		}

		//calc  temp_vector_transpose[] *  Sk_transpose[]
		//PASS
		//cout << "\nTOP MATRIX VALUE:\n";
		for (int row = 0; row < 2; row++) {
			
			for (int col = 0; col < 2; col++) {
				temp_matrix[row][col] = 0;
				temp_matrix[row][col] = (temp_vector_transpose[row] * Sk_transpose[col]);
				//cout << temp_matrix[row][col] << "\t";
			}
			//cout << "\n";
		}

		//multiply top by 1/bottom_Val
		//PASS
		//cout << "\n1/2 * TOP MATRIX VALUE:\n";
		for (int row = 0; row < 2; row++) {
			for (int col = 0; col < 2; col++) {

				//temp_matrix[row][col] = ((1/bottom_Val) * temp_matrix[row][col]);
				temp_matrix[row][col] *= (1 / bottom_Val);
				//cout << temp_matrix[row][col] << "\t";
			}
			//cout << "\n";
		}

		//PASS
		//cout << "\nNEW JACOBIAN:\n";

		//final step Bk + the resulting matrix calculated above
		for (int row = 0; row < 2; row++) {
			for (int col = 0; col < 2; col++) {

				//new Jacobian matrix or Bk+1
				jacobian[row][col] = (jacobian[row][col] + temp_matrix[row][col]);

				//cout << jacobian[row][col] << "\t";

			}
			//cout << "\n";
		}

		//break;


		//Final Step
		//set current Xk to new Xk before starting next iteration of loop
		for (int row = 0; row < 2; row++) {

			Xk_transpose[row] = new_Xk_transpose[row];

		}

		k++;
		//end of algorithm

		if ((calcTwoNorm(temp_NEW_Fx_vector_result, 2) - calcTwoNorm(temp_Fx_vector_result, 2))/ calcTwoNorm(temp_NEW_Fx_vector_result, 2)*100 <= 0.05) {
			return 0;
		}

	}

	return 0;
}