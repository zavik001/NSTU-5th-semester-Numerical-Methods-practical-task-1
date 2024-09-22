#include "SLAUSolverLDLT.hpp"

SLAUSolverLDLT::SLAUSolverLDLT(const string &inputFilePath,
                               const string &alFilePath,
                               const string &dFilePath,
                               const string &fFilePath,
                               const string &solveFilePath)
    : solveFilePath(solveFilePath), AlFilePath(alFilePath), DFilePath(dFilePath)
{
    cout << fixed << setprecision(PRECISION_DIGITS);

    loadFromFile(inputFilePath);

    matrixAL.resize(n, vector<floatingPointType>(m, 0.0));
    diagD.resize(n, 0.0);
    vectorF.resize(n, 0.0);

    loadFromFile(alFilePath, matrixAL);
    loadFromFile(dFilePath, diagD);
    loadFromFile(fFilePath, vectorF);
}

void SLAUSolverLDLT::loadFromFile(const string &filePath)
{
    ifstream file(filePath);
    if (!file.is_open())
    {
        throw runtime_error("Could not open file: " + filePath);
    }
    file >> n >> m;
    file.close();
}

void SLAUSolverLDLT::loadFromFile(const string &filePath, vector<vector<floatingPointType>> &matrix)
{
    ifstream file(filePath);
    if (!file.is_open())
    {
        throw runtime_error("Could not open file: " + filePath);
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            file >> matrix[i][j];
        }
    }
    file.close();
}

void SLAUSolverLDLT::loadFromFile(const string &filePath, vector<floatingPointType> &vector)
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

void SLAUSolverLDLT::performLDLtDecomposition()
{
    for (int i = 0; i < n; ++i)
    {
        sum sumD = 0.0;

        for (int j = max(0, i - m); j <= i + m && j < n; ++j)
        {
            if (j < i)
            {
                sumD += matrixAL[i][m - i + j] * matrixAL[i][m - i + j] * diagD[j];
            }
            else if (j == i)
            {
                diagD[i] -= sumD;
            }
            else if (j > i)
            {
                sum sumL = 0.0;
                for (int g = max(0, i - m); g < i; ++g)
                {
                    sumL += matrixAL[j][m - j + g] * matrixAL[i][m - i + g] * diagD[g];
                }
                matrixAL[j][m - j + i] = (matrixAL[j][m - j + i] - sumL) / diagD[i];
            }
        }
    }
}

void SLAUSolverLDLT::solveForwardSubstitution()
{
    for (int i = 0; i < n; ++i)
    {
        sum sum = 0.0;
        for (int j = max(0, i - m); j < i; ++j)
        {
            sum += matrixAL[i][m - i + j] * vectorF[j];
        }
        vectorF[i] -= sum;
    }
}

void SLAUSolverLDLT::solveDiagonalSubstitution()
{
    for (int i = 0; i < n; ++i)
    {
        vectorF[i] /= diagD[i];
    }
}

void SLAUSolverLDLT::solveBackwardSubstitution()
{
    for (int i = n - 1; i >= 0; --i)
    {
        floatingPointType sum = 0.0;
        for (int j = i + 1; j < n && j <= i + m; ++j)
        {
            sum += matrixAL[j][m - j + i] * vectorF[j];
        }
        vectorF[i] -= sum;
    }
}

void SLAUSolverLDLT::solveLinearSystem()
{
    solveForwardSubstitution();
    solveDiagonalSubstitution();
    solveBackwardSubstitution();
}

void SLAUSolverLDLT::writeSolutionToFile()
{
    ofstream outFile(solveFilePath);
    if (!outFile.is_open())
    {
        throw runtime_error("Could not open file: " + solveFilePath);
    }

    outFile << fixed << setprecision(PRECISION_DIGITS);
    for (const auto &value : vectorF)
    {
        outFile << value << '\n';
    }

    outFile.close();
}

void SLAUSolverLDLT::printVector()
{
    for (int i = 0; i < n; i++)
    {
        cout << vectorF[i] << '\n';
    }
    cout << '\n';
}

void SLAUSolverLDLT::returnMatix()
{
    loadFromFile(AlFilePath, matrixAL);
    loadFromFile(DFilePath, diagD);
}

void SLAUSolverLDLT::printMultiplyMatrixToVector()
{
    cout << "Result of multiplying matrix (A = AL + D) by vector X:" << "\n";

    for (int i = 0; i < n; ++i)
    {
        floatingPointType result = diagD[i] * vectorF[i];

        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                continue;
            }
            else if (j < i && (i - j) <= m)
            {
                result += matrixAL[i][m - i + j] * vectorF[j];
            }
            else if (i < j && (j - i) <= m)
            {
                result += matrixAL[j][m - j + i] * vectorF[j];
            }
        }

        cout << result << ' ';
    }

    cout << endl;
}

void SLAUSolverLDLT::printRestoredMatrix()
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                cout << diagD[i] << " ";
            }
            else if (j < i && (i - j) <= m)
            {
                cout << matrixAL[i][m - i + j] << " ";
            }
            else if (i < j && (j - i) <= m)
            {
                cout << matrixAL[j][m - j + i] << " ";
            }
            else
            {
                cout << 0.0 << " ";
            }
        }
        cout << '\n';
    }
    cout << '\n';
}

void SLAUSolverLDLT::printMatrixAL()
{
    for (const auto &row : matrixAL)
    {
        for (const auto &value : row)
        {
            cout << value << " ";
        }
        cout << '\n';
    }
    cout << '\n';
}