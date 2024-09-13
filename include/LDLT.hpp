#ifndef LDLT_HPP
#define LDLT_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

using namespace std;

#if defined(FF)
typedef float floatingPointType;
#elif defined(DD)
typedef double floatingPointType;
#elif defined(LDD)
typedef long double floatingPointType;
#endif

class LDLT
{
private:
    vector<vector<floatingPointType>> matrixAL; // Matrix for decomposition
    vector<floatingPointType> diagD;            // Diagonal elements
    vector<floatingPointType> vectorF;          // Right-hand side vector
    vector<floatingPointType> y, z, X, result;  // Temporary vectors for solving the system

    int n; // Matrix size
    int k; // Bandwidth

    string solveFilePath; // Path to save solution (x.txt)

public:
    // Constructor to initialize dimensions, reserve memory and load data
    LDLT(const string &inputFilePath,
         const string &alFilePath,
         const string &dFilePath,
         const string &fFilePath,
         const string &outputFilePath);

    // Load matrix data from file
    void loadFromFile(const string &filePath, vector<vector<floatingPointType>> &matrix);

    // Load vector data from file
    void loadFromFile(const string &filePath, vector<floatingPointType> &vector);

    // Load matrix dimensions (n and k) from file
    void loadFromFile(const string &filePath);

    // Print matrix AL (lower triangular matrix)
    void printMatrixAL();

    // Print the restored matrix n*k + d => n*n
    void printRestoredMatrix();

    // Perform LDL^T decomposition on the matrix
    void performLDLtDecomposition();

    // Multiply band matrix by itself to restore the matrix A from L, D, and L^T
    void multiplyBandMatrix();

    // Solve the system of linear equations
    void solveLinearSystem();

    // Print result of matrix-vector multiplication
    void printMultiplyByVector();

    // Write the solution (vector X) to a file
    void writeSolutionToFile();
};

#endif // LDLT_HPP