#include "LDLT.hpp"

// Constructor that initializes matrix, diagonal, vector F, and reserve memory for solution
LDLT::LDLT(const string &inputFilePath,
           const string &alFilePath,
           const string &dFilePath,
           const string &fFilePath,
           const string &solveFilePath)
    : solveFilePath(solveFilePath), AlFilePath(alFilePath), DFilePath(dFilePath)
{
    // Load matrix dimensions
    loadFromFile(inputFilePath);

    // Resize vectors and matrix for the problem dimensions
    matrixAL.resize(n, vector<floatingPointType>(k, 0.0));
    diagD.resize(n, 0.0);
    vectorF.resize(n, 0.0);
    vectorX.resize(n, 0.0); // Solution vector
    y.resize(n, 0.0); // Intermediate vector y
    z.resize(n, 0.0); // Intermediate vector z

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
            sum += matrixAL[j][k - j + i] * vectorX[j];
        }
        vectorX[i] = z[i] - sum;
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

// Print result of matrix-vector multiplication
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
