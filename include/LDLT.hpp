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

/**
 * @class LDLT
 * @brief This class implements the LDL^T decomposition algorithm for 
 * symmetric matrices, specifically for indefinite (non-positive definite) matrices.
 * The matrix A is decomposed as A = L * D * L^T, where:
 * - L is a lower triangular matrix,
 * - D is a diagonal matrix with ±1 on the diagonal.
 */
class LDLT
{
private:
    vector<vector<floatingPointType>> matrixAL; // Matrix for decomposition (lower triangular part)
    vector<floatingPointType> diagD;            // Diagonal elements (D)
    vector<floatingPointType> vectorF;          // Right-hand side vector (F)
    vector<floatingPointType> vectorX;          // Solution vector (X)
    vector<floatingPointType> y, z;             // Intermediate vectors for solving the system

    int n; // Size of the matrix
    int k; // Bandwidth of the matrix

    string solveFilePath; // Path to save the solution (x.txt)
    string AlFilePath;    // Path to the file containing matrix A
    string DFilePath;     // Path to the file containing diagonal D

public:
    /**
     * @brief Constructor that initializes the dimensions of the matrix,
     * allocates memory, and loads the data from files.
     * 
     * @param inputFilePath File path containing the matrix dimensions.
     * @param alFilePath File path for the matrix AL.
     * @param dFilePath File path for the diagonal matrix D.
     * @param fFilePath File path for the right-hand side vector F.
     * @param outputFilePath File path to save the solution.
     */
    LDLT(const string &inputFilePath,
         const string &alFilePath,
         const string &dFilePath,
         const string &fFilePath,
         const string &outputFilePath);

    /**
     * @brief Loads matrix dimensions (n and k) from a file.
     * 
     * @param filePath Path to the file with matrix dimensions.
     */
    void loadFromFile(const string &filePath);

    /**
     * @brief Loads matrix data from a file into a 2D vector.
     * 
     * @param filePath Path to the file with matrix data.
     * @param matrix Matrix to store the loaded data.
     */
    void loadFromFile(const string &filePath, vector<vector<floatingPointType>> &matrix);

    /**
     * @brief Loads vector data from a file into a 1D vector.
     * 
     * @param filePath Path to the file with vector data.
     * @param vector Vector to store the loaded data.
     */
    void loadFromFile(const string &filePath, vector<floatingPointType> &vector);

    /**
     * @brief Performs LDL^T decomposition of matrix A.
     * Decomposes the symmetric matrix A into L (lower triangular), 
     * D (diagonal matrix), and L^T (transpose of L).
     * 
     * Time complexity: O(n * k), where n is the size of the matrix and 
     * k is the bandwidth.
     */
    void performLDLtDecomposition();

    /**
     * @brief Solves the system of linear equations using LDL^T decomposition.
     * The solution is obtained in three phases:
     * 1. Forward substitution for Ly = F.
     * 2. Diagonal substitution for Dz = y.
     * 3. Backward substitution for L^T X = z.
     * 
     * Time complexity: O(n * k), where n is the size of the matrix and 
     * k is the bandwidth.
     */
    void solveLinearSystem();

    /**
     * @brief Writes the solution vector X to a file.
     */
    void writeSolutionToFile();

    /**
     * @brief Restores the original matrix A by multiplying L, D, and L^T.
     * Multiplies the band matrix by itself to recover matrix A.
     */
    void returnMatix();

    /**
     * @brief Prints the result of multiplying matrix A by vector X.
     */
    void printMultiplyMatrixToVectorX();

    /**
     * @brief Prints the lower triangular matrix AL.
     */
    void printMatrixAL();

    /**
     * @brief Prints the restored matrix A (after decomposition).
     */
    void printRestoredMatrix();

    /**
     * @brief Prints vectors X (solution) and F (right-hand side vector).
     */
    void printVectors();
};

#endif // LDLT_HPP
