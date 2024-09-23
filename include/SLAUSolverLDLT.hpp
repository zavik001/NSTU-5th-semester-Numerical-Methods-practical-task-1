/**
 * @file SLAUSolverLDLT.hpp
 * @brief Header file for solving systems of linear algebraic equations (SLAE) using the LDLT decomposition.
 *
 * This class implements the LDLT decomposition for symmetric matrices stored in a banded format,
 * and methods for solving the resulting linear systems.
 */

#ifndef SLAUSolverLDLT_HPP
#define SLAUSolverLDLT_HPP

#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

#if defined(FF)
typedef float floatingPointType;
typedef float sum;
constexpr int PRECISION_DIGITS = 7;
#elif defined(DD)
typedef double floatingPointType;
typedef double sum;
constexpr int PRECISION_DIGITS = 15;
#elif defined(FD)
typedef float floatingPointType;
typedef double sum;
constexpr int PRECISION_DIGITS = 7;
#endif

/**
 * @class SLAUSolverLDLT
 * @brief A class for solving SLAE using LDLT decomposition with a banded matrix format.
 */
class SLAUSolverLDLT
{
private:
    vector<vector<floatingPointType>> matrixAL; ///< The lower triangular matrix in banded form (L)
    vector<floatingPointType> diagD;            ///< The diagonal matrix (D)
    vector<floatingPointType> vectorF;          ///< Vector used for solving the system

    int n; ///< The size of the system (number of equations)
    int m; ///< The bandwidth of the matrix

    string solveFilePath; ///< Path to the file for output results
    string AlFilePath;    ///< Path to the file containing matrix A (banded part)
    string DFilePath;     ///< Path to the file containing diagonal matrix D

public:
    /**
     * @brief Initializes the matrix size and allocates memory.
     *
     * @param a Number of equations (matrix size)
     * @param b Matrix bandwidth
     */
    void initialize(int a, int b);

    /**
     * @brief Constructor that loads matrix and vector data from files.
     *
     * @param inputFilePath Path to the input file with system size
     * @param alFilePath Path to the file with matrix A (banded part)
     * @param dFilePath Path to the file with matrix D
     * @param fFilePath Path to the file with the vector F
     * @param outputFilePath Path to the file for saving the solution
     */
    SLAUSolverLDLT(const string &inputFilePath,
                   const string &alFilePath,
                   const string &dFilePath,
                   const string &fFilePath,
                   const string &outputFilePath);

    /**
     * @brief Loads matrix dimensions from a file.
     *
     * @param filePath Path to the file with matrix dimensions
     * @param a Reference to store the number of rows
     * @param b Reference to store the matrix bandwidth
     */
    void loadFromFile(const string &filePath, int &a, int &b);

    /**
     * @brief Loads a matrix from a file.
     *
     * @param filePath Path to the file containing the matrix
     * @param matrix Reference to the matrix where data will be stored
     */
    void loadFromFile(const string &filePath, vector<vector<floatingPointType>> &matrix);

    /**
     * @brief Loads a vector from a file.
     *
     * @param filePath Path to the file containing the vector
     * @param vector Reference to the vector where data will be stored
     */
    void loadFromFile(const string &filePath, vector<floatingPointType> &vector);

    /**
     * @brief Performs LDLT decomposition of the matrix.
     *
     * Decomposes the matrix A into L, D, and L^T where:
     * L - Lower triangular matrix
     * D - Diagonal matrix A
     */
    void performLDLtDecomposition();

    /**
     * @brief Solves the system using forward substitution for L * y = b.
     *
     * This method solves the first step of the LDLT process, where the matrix L is
     * the lower triangular part and y is an intermediate solution vector.
     */
    void solveForwardSubstitution();

    /**
     * @brief Solves the system using diagonal substitution for D * z = y.
     *
     * This method solves the second step of the LDLT process, where the matrix D
     * is the diagonal matrix, and z is the solution vector after the diagonal substitution.
     */
    void solveDiagonalSubstitution();

    /**
     * @brief Solves the system using backward substitution for L^T * x = z.
     *
     * This method solves the final step of the LDLT process, where L^T is the transpose
     * of the lower triangular matrix, and x is the final solution vector for the system.
     */
    void solveBackwardSubstitution();

    /**
     * @brief Solves the full linear system using LDLT decomposition.
     *
     * This method combines forward substitution, diagonal substitution, and backward
     * substitution to find the solution of the system of linear equations.
     */
    void solveLinearSystem();

    /**
     * @brief Writes the solution vector F to a file.
     *
     * Saves the solution of the system to the specified file.
     */
    void writeVectorFToFile();

    /**
     * @brief Prints the vector F to the console.
     *
     * Outputs the current state of the vector F, which contains the solution.
     */
    void printvectorF();

    /**
     * @brief Reloads the matrix and diagonal values from files.
     *
     * Re-initializes matrixAL and diagD by loading their values from files.
     */
    void returnMatix();

    /**
     * @brief Multiplies the restored matrix (A = L + D) by the solution vector and prints the result.
     *
     * Computes the product of matrix A and vector X, where matrix A is reconstructed from
     * its banded format representation, and prints the resulting vector.
     */
    void printMultiplyMatrixToVector();

    /**
     * @brief Prints the fully restored matrix (A = L + D) to the console.
     *
     * Reconstructs and prints the original matrix A, which is the sum of the lower triangular
     * matrix L and the diagonal matrix D, from the banded format.
     */
    void printRestoredMatrix();

    /**
     * @brief Prints the lower triangular matrix (L) stored in banded format.
     *
     * Outputs the matrixAL to the console, showing the banded lower triangular part of the matrix.
     */
    void printMatrixAL();

    /**
     * @brief Generates a Hilbert matrix in banded format and initializes vector F.
     *
     * Fills the matrixAL with values from a Hilbert matrix (used for testing numerical stability)
     * and initializes the vectorF with values.
     */
    void HilbertBandMatrix();

    /**
     * @brief Restores the Hilbert matrix in banded format after modifications.
     *
     * Recalculates the matrixAL and diagD using values from a Hilbert matrix.
     */
    void returnMatixAfterHilbert();
};

#endif // SLAUSolverLDLT_HPP
