#include "LDLT.hpp"

LDLT::LDLT(const string &inputFilePath,
           const string &alFilePath,
           const string &dFilePath,
           const string &fFilePath,
           const string &solveFilePath)
    : solveFilePath(solveFilePath), AlFilePath(alFilePath), DFilePath(dFilePath)
{
    loadFromFile(inputFilePath);

    matrixAL.resize(n, vector<floatingPointType>(k, 0.0));
    diagD.resize(n, 0.0);
    vectorF.resize(n, 0.0);
    vectorX.resize(n, 0.0);
    y.resize(n, 0.0);
    z.resize(n, 0.0);

    loadFromFile(alFilePath, matrixAL);
    loadFromFile(dFilePath, diagD);
    loadFromFile(fFilePath, vectorF);
}

void LDLT::loadFromFile(const string &filePath)
{
    ifstream file(filePath);
    if (!file.is_open())
    {
        throw runtime_error("Could not open file: " + filePath);
    }
    file >> n >> k;
    file.close();
}

void LDLT::loadFromFile(const string &filePath, vector<vector<floatingPointType>> &matrix)
{
    ifstream file(filePath);
    if (!file.is_open())
    {
        throw runtime_error("Could not open file: " + filePath);
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < k; ++j)
        {
            file >> matrix[i][j];
        }
    }
    file.close();
}

void LDLT::loadFromFile(const string &filePath, vector<floatingPointType> &vector)
{
    ifstream file(filePath);
    if (!file.is_open())
    {
        throw runtime_error("Could not open file: " + filePath);
    }
    for (int i = 0; i < n; ++i)
    {
        file >> vector[i];
    }
    file.close();
}

/**
 * @brief Perform LDL^T decomposition on the matrix A.
 *
 * This method performs the decomposition of the matrix A into
 * a lower triangular matrix (L), a diagonal matrix (D), and the transpose of L (L^T).
 * The algorithm has a time complexity of O(nk) due to the band structure of the matrix.
 * OpenMP is used to parallelize the outer loop to improve performance, especially for large matrices.
 *
 * @note The decomposition modifies the matrix A in place.
 */
void LDLT::performLDLtDecomposition()
{
// Parallelizing the outer loop with OpenMP
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        floatingPointType sumD = 0.0;

        // Compute the sum for D[i] using previously computed values in L and D
        for (int j = max(0, i - k); j < i; ++j)
        {
            sumD += matrixAL[i][k - i + j] * matrixAL[i][k - i + j] * diagD[j];
        }
        diagD[i] -= sumD;

        // Update the L[j][i] elements for subsequent rows
        for (int j = i + 1; j < n && j <= i + k; ++j)
        {
            floatingPointType sumL = 0.0;
            for (int m = max(0, i - k); m < i; ++m)
            {
                sumL += matrixAL[j][k - j + m] * matrixAL[i][k - i + m] * diagD[m];
            }
            matrixAL[j][k - j + i] = (matrixAL[j][k - j + i] - sumL) / diagD[i];
        }
    }
}

/**
 * @brief Solve the system of linear equations using the LDL^T decomposition.
 *
 * This method uses forward substitution, diagonal substitution,
 * and backward substitution to solve the system of equations A * X = F.
 *
 * Forward substitution and backward substitution are inherently sequential due to data dependencies,
 * but diagonal substitution is parallelized with OpenMP, as each element can be computed independently.
 *
 * Time complexity is O(nk) due to the band structure of the matrix.
 */
void LDLT::solveLinearSystem()
{
    // Forward substitution: L * y = F
    // This loop cannot be parallelized due to the sequential nature of forward substitution.
    for (int i = 0; i < n; ++i)
    {
        floatingPointType sum = 0.0;
        for (int j = max(0, i - k); j < i; ++j)
        {
            sum += matrixAL[i][k - i + j] * y[j];
        }
        y[i] = vectorF[i] - sum;
    }

// Diagonal substitution: D * z = y
// This loop can be parallelized as each element is independent.
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        z[i] = y[i] / diagD[i];
    }

    // Backward substitution: L^T * X = z
    // This loop also cannot be parallelized due to sequential data dependencies.
    for (int i = n - 1; i >= 0; --i)
    {
        floatingPointType sum = 0.0;
        for (int j = i + 1; j < n && j <= i + k; ++j)
        {
            sum += matrixAL[j][k - j + i] * vectorX[j];
        }
        vectorX[i] = z[i] - sum;
    }

    // Save the solution to a file
    writeSolutionToFile();
}

void LDLT::printMatrixAL()
{
    for (const auto &row : matrixAL)
    {
        for (const auto &value : row)
        {
            cout << setw(10) << value << " ";
        }
        cout << '\n';
    }
    cout << '\n';
}

void LDLT::printRestoredMatrix()
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                cout << setw(10) << diagD[i] << " ";
            }
            else if (j < i && (i - j) <= k)
            {
                cout << setw(10) << matrixAL[i][k - i + j] << " ";
            }
            else if (i < j && (j - i) <= k)
            {
                cout << setw(10) << matrixAL[j][k - j + i] << " ";
            }
            else
            {
                cout << setw(10) << 0.0 << " ";
            }
        }
        cout << '\n';
    }
    cout << '\n';
}

void LDLT::writeSolutionToFile()
{
    ofstream outFile(solveFilePath);
    if (!outFile.is_open())
    {
        throw runtime_error("Could not open file: " + solveFilePath);
    }

    for (const auto &value : vectorX)
    {
        outFile << value << '\n';
    }

    outFile.close();
}

void LDLT::printVectors()
{
    cout << "X:\n";
    for (int i = 0; i < n; i++)
    {
        cout << vectorX[i] << '\n';
    }
    cout << '\n';

    cout << "F:\n";
    for (int i = 0; i < n; i++)
    {
        cout << vectorF[i] << '\n';
    }
    cout << '\n';
}

void LDLT::returnMatix()
{
    loadFromFile(AlFilePath, matrixAL);
    loadFromFile(DFilePath, diagD);
}

void LDLT::printMultiplyMatrixToVectorX()
{
    vector<floatingPointType> result(n, 0.0);

    for (int i = 0; i < n; ++i)
    {
        result[i] += diagD[i] * vectorX[i];

        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                continue;
            }
            else if (j < i && (i - j) <= k)
            {
                result[i] += matrixAL[i][k - i + j] * vectorX[j];
            }
            else if (i < j && (j - i) <= k)
            {
                result[i] += matrixAL[j][k - j + i] * vectorX[j];
            }
        }
    }

    cout << "Result of multiplying matrix (A = AL + D) by vector X:" << endl;
    for (int i = 0; i < n; ++i)
    {
        cout << setprecision(7) << result[i] << ' ';
    }
    cout << endl;
}
