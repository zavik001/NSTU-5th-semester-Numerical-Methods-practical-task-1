#include "LDLT.hpp"

// Constructor that initializes matrix, diagonal, vector F, and reserve memory for solution
LDLT::LDLT(const string &inputFilePath,
           const string &alFilePath,
           const string &dFilePath,
           const string &fFilePath,
           const string &solveFilePath)
    : solveFilePath(solveFilePath)
{
    // Load matrix dimensions
    loadFromFile(inputFilePath);

    // Resize vectors and matrix for the problem dimensions
    matrixAL.resize(n, vector<floatingPointType>(k, 0.0));
    diagD.resize(n, 0.0);
    vectorF.resize(n, 0.0);
    X.resize(n, 0.0); // Solution vector
    y.resize(n, 0.0); // Intermediate vector y
    z.resize(n, 0.0); // Intermediate vector z
    result.resize(n, 0.0);

    // Load matrix AL, diagonal D, and vector F from the specified files
    loadFromFile(alFilePath, matrixAL);
    loadFromFile(dFilePath, diagD);
    loadFromFile(fFilePath, vectorF);
}

// Load matrix dimensions from a file (n and k)
void LDLT::loadFromFile(const string &filePath)
{
    ifstream file(filePath);
    if (!file.is_open())
    {
        throw runtime_error("Could not open file: " + filePath);
    }
    file >> n >> k; // Read matrix dimensions
    file.close();
}

// Load matrix data from a file into a vector (matrix)
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

// Load vector data from a file into a 1D vector
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

// Print matrix AL in a formatted way
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

// Print restored matrix after LDL^T decomposition
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

// Perform LDL^T decomposition
void LDLT::performLDLtDecomposition()
{
    for (int i = 0; i < n; ++i)
    {
        floatingPointType sumD = 0.0;
        for (int j = max(0, i - k); j < i; ++j)
        {
            sumD += matrixAL[i][k - i + j] * matrixAL[i][k - i + j] * diagD[j];
        }
        diagD[i] = diagD[i] - sumD;

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

void LDLT::multiplyBandMatrix()
{
    // Step 1: Create matrix A to store L * D * L^T
    vector<vector<floatingPointType>> A(n, vector<floatingPointType>(n, 0.0));

    // Step 2: Compute L * D * L^T
    for (int i = 0; i < n; i++)
    {
        // Diagonal elements of A
        A[i][i] = diagD[i]; // Only diagonal D[i]

        // Below-diagonal elements
        for (int j = max(0, i - k); j < i; ++j)
        {
            floatingPointType sum = 0.0;
            for (int m = max(0, j - k); m < j; ++m)
            {
                sum += matrixAL[i][k - i + m] * diagD[m] * matrixAL[j][k - j + m]; // L[i,m] * D[m] * L[j,m]
            }
            A[i][j] = matrixAL[i][k - i + j] * diagD[j] + sum; // L[i,j] * D[j]
            A[j][i] = A[i][j];                                 // Symmetry for upper triangular
        }
    }

    // Step 3: Output restored matrix A
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << setw(10) << A[i][j] << " ";
        }
        cout << '\n';
    }
}

// Solve the system using LDL^T decomposition
void LDLT::solveLinearSystem()
{
    // Forward substitution for Ly = F
    for (int i = 0; i < n; ++i)
    {
        floatingPointType sum = 0.0;
        for (int j = max(0, i - k); j < i; ++j)
        {
            sum += matrixAL[i][k - i + j] * y[j];
        }
        y[i] = (vectorF[i] - sum);
    }

    // Diagonal substitution for Dz = y
    for (int i = 0; i < n; ++i)
    {
        z[i] = y[i] / diagD[i];
    }

    // Backward substitution for L^T X = z
    for (int i = n - 1; i >= 0; --i)
    {
        floatingPointType sum = 0.0;
        for (int j = i + 1; j < n && j <= i + k; ++j)
        {
            sum += matrixAL[j][k - j + i] * X[j];
        }
        X[i] = z[i] - sum;
    }

    // Write solution to file
    writeSolutionToFile();
}

// Write solution to a file (vector X)
void LDLT::writeSolutionToFile()
{
    ofstream outFile(solveFilePath);
    if (!outFile.is_open())
    {
        throw runtime_error("Could not open file: " + solveFilePath);
    }

    for (const auto &value : X)
    {
        outFile << value << '\n';
    }

    outFile.close();
}

// Print result of matrix-vector multiplication
void LDLT::printMultiplyByVector()
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = max(0, i - k); j <= i; ++j)
        {
            result[i] += matrixAL[i][k - i + j] * X[j];
            if (i != j)
            {
                result[j] += matrixAL[i][k - i + j] * X[i];
            }
        }
    }

    for (int i = 0; i < n; ++i)
    {
        result[i] += diagD[i] * X[i];
    }

    for (const auto &value : result)
    {
        cout << value << '\n';
    }
}