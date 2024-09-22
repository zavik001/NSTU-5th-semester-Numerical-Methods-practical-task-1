/**
 * @file SLAUSolverLDLT.hpp
 * @brief Header file for the SLAUSolverLDLT class.
 *
 * Contains the declaration of the LDLT decomposition solver with banded matrix storage.
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
 * @brief Class for solving systems of linear equations (SLAE) using LDLT decomposition.
 *
 * The matrix is stored in banded format. This class performs LDLT decomposition,
 * where L is a lower triangular matrix with 1's on the diagonal, and D is a diagonal matrix.
 */
class SLAUSolverLDLT
{
private:
    /// Banded format matrix A (stored in L) as a 2D vector
    vector<vector<floatingPointType>> matrixAL;

    /// Diagonal elements of matrix D
    vector<floatingPointType> diagD;

    /// Vector F representing the right-hand side of the equation Ax = b
    vector<floatingPointType> vectorF;

    /// Number of equations (n) and bandwidth (m)
    int n;
    int m;

    /// File paths for storing results and loading matrices
    string solveFilePath;
    string AlFilePath;
    string DFilePath;

public:
    /**
     * @brief Constructor that initializes the solver with input file paths.
     *
     * @param inputFilePath Path to the file containing matrix dimensions.
     * @param alFilePath Path to the file containing the lower triangular matrix (banded format).
     * @param dFilePath Path to the file containing the diagonal matrix.
     * @param fFilePath Path to the file containing the right-hand side vector.
     * @param outputFilePath Path to the file where the solution will be written.
     */
    SLAUSolverLDLT(const string &inputFilePath,
                   const string &alFilePath,
                   const string &dFilePath,
                   const string &fFilePath,
                   const string &outputFilePath);

    /**
     * @brief Loads matrix or vector data from a file.
     *
     * This method reads the matrix dimensions from the file.
     * @param filePath The path to the file containing the matrix/vector data.
     */
    void loadFromFile(const string &filePath);

    /**
     * @brief Loads matrix data from a file into a 2D vector.
     *
     * @param filePath Path to the file containing the matrix data.
     * @param matrix 2D vector where the matrix data will be stored.
     */
    void loadFromFile(const string &filePath, vector<vector<floatingPointType>> &matrix);

    /**
     * @brief Loads vector data from a file into a 1D vector.
     *
     * @param filePath Path to the file containing the vector data.
     * @param vector 1D vector where the data will be stored.
     */
    void loadFromFile(const string &filePath, vector<floatingPointType> &vector);

    /**
     * @brief Performs the LDLT decomposition of the matrix A.
     *
     * The method computes the lower triangular matrix L and diagonal matrix D.
     */
    void performLDLtDecomposition();

    /**
     * @brief Solves the system using forward substitution.
     *
     * This method solves the system Ly = b, where L is lower triangular.
     */
    void solveForwardSubstitution();

    /**
     * @brief Solves the system with the diagonal matrix D.
     *
     * This method solves the system Dz = y, where D is diagonal.
     */
    void solveDiagonalSubstitution();

    /**
     * @brief Solves the system using backward substitution.
     *
     * This method solves the system L^T x = z, where L^T is the transpose of L.
     */
    void solveBackwardSubstitution();

    /**
     * @brief Solves the complete linear system Ax = b.
     *
     * This method calls forward, diagonal, and backward substitutions in sequence.
     */
    void solveLinearSystem();

    /**
     * @brief Writes the solution vector x to a file.
     *
     * The solution is written to the file specified by solveFilePath.
     */
    void writeSolutionToFile();

    /**
     * @brief Restores the matrix from files.
     *
     * This method reloads matrix A (stored in L and D) from the corresponding files.
     */
    void returnMatix();

    /**
     * @brief Prints the result of multiplying the restored matrix A by the solution vector X.
     */
    void printMultiplyMatrixToVector();

    /**
     * @brief Prints the solution vector.
     */
    void printVector();

    /**
     * @brief Prints the banded matrix AL.
     */
    void printMatrixAL();

    /**
     * @brief Prints the restored matrix A from L and D.
     */
    void printRestoredMatrix();
};

#endif // SLAUSolverLDLT_HPP
