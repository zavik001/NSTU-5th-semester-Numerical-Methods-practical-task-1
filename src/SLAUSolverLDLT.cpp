/**
 * @file SLAUSolverLDLT.cpp
 * @brief Implementation of the SLAUSolverLDLT class for solving 
 *        systems of linear algebraic equations (SLAE) using LDL^T decomposition.
 *
 * This file contains the implementation of methods for the SLAUSolverLDLT class,
 * which uses the LDL^T decomposition to solve symmetric banded matrices.
 * 
 * The main functionalities include:
 * - LDL^T decomposition of the matrix
 * - Solving systems of linear equations with forward, diagonal, and backward substitution
 * - Matrix and vector I/O operations
 * - Printing the results for validation and debugging
 * 
 * The class supports banded matrix storage for optimized memory usage and performance.
 *
 * @author zavik001
 * @version 1.0
 * 
 * Usage example:
 * @code
 * SLAUSolverLDLT solver("input.txt", "al.txt", "d.txt", "f.txt", "output.txt");
 * solver.performLDLtDecomposition();
 * solver.solveLinearSystem();
 * solver.writeVectorFToFile();
 * @endcode
 */
#include "SLAUSolverLDLT.hpp"

void SLAUSolverLDLT::initialize(int a, int b)
{
    n = a;
    m = b;
    matrixAL.resize(n, vector<floatingPointType>(m, 0.0));
    diagD.resize(n, 0.0);
    vectorF.resize(n, 0.0);
}

SLAUSolverLDLT::SLAUSolverLDLT(const string &inputFilePath,
                               const string &alFilePath,
                               const string &dFilePath,
                               const string &fFilePath,
                               const string &solveFilePath)
    : solveFilePath(solveFilePath), AlFilePath(alFilePath), DFilePath(dFilePath)
{
    cout << fixed << setprecision(PRECISION_DIGITS);

    int a = 0, b = 0;
    loadFromFile(inputFilePath, a, b);
    initialize(a, b);

    loadFromFile(alFilePath, matrixAL);
    loadFromFile(dFilePath, diagD);
    loadFromFile(fFilePath, vectorF);
}

void SLAUSolverLDLT::loadFromFile(const string &filePath, int &a, int &b)
{
    ifstream file(filePath);
    if (!file.is_open())
    {
        throw runtime_error("Could not open file: " + filePath);
    }
    file >> a >> b;
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
        floatingPointType sumD = 0;
        int baseIndexI = m - i;

        for (int j = max(0, i - m); j <= i + m && j < n; ++j)
        {
            int baseIndexJ = m - j;

            if (j < i)
            {
                int indexIJ = baseIndexI + j;
                sumD += matrixAL[i][indexIJ] * matrixAL[i][indexIJ] * diagD[j];
            }
            else if (j == i)
            {
                diagD[i] -= sumD;
            }
            else if (j > i)
            {
                floatingPointType sumL = 0;

                int indexIK = baseIndexI;
                int indexJK = baseIndexJ;
                for (int k = max(0, i - m); k < i; ++k)
                {
                    indexIK += k;
                    indexJK += k;
                    sumL += matrixAL[j][indexJK] * matrixAL[i][indexIK] * diagD[k];
                }

                int indexJI = baseIndexJ + i;
                matrixAL[j][indexJI] = (matrixAL[j][indexJI] - sumL) / diagD[i];
            }
        }
    }
}

void SLAUSolverLDLT::solveForwardSubstitution()
{
    for (int i = 0; i < n; ++i)
    {
        floatingPointType sum = 0;
        int baseIndexI = m - i;

        for (int j = max(0, i - m); j < i; ++j)
        {
            int indexIJ = baseIndexI + j;
            sum += matrixAL[i][indexIJ] * vectorF[j];
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
        int baseIndexI = m - i;

        for (int j = i + 1; j < n && j <= i + m; ++j)
        {
            int baseIndexJ = m - j;
            int indexJI = baseIndexJ + i;
            sum += matrixAL[j][indexJI] * vectorF[j];
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

void SLAUSolverLDLT::writeVectorFToFile()
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

void SLAUSolverLDLT::printvectorF()
{
    for (const auto &val : vectorF)
    {
        cout << val << '\n';
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
        int baseIndexI = m - i;

        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                continue;
            }
            else if (j < i && (i - j) <= m)
            {
                int indexIJ = baseIndexI + j;
                result += matrixAL[i][indexIJ] * vectorF[j];
            }
            else if (i < j && (j - i) <= m)
            {
                int baseIndexJ = m - j;
                int indexJI = baseIndexJ + i;
                result += matrixAL[j][indexJI] * vectorF[j];
            }
        }
        cout << result << endl;
    }
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
                int baseIndexI = m - i;
                int indexIJ = baseIndexI + j;
                cout << matrixAL[i][indexIJ] << " ";
            }
            else if (i < j && (j - i) <= m)
            {
                int baseIndexJ = m - j;
                int indexJI = baseIndexJ + i;
                cout << matrixAL[j][indexJI] << " ";
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

void SLAUSolverLDLT::HilbertBandMatrix()
{
    for (int i = 1; i < n; ++i)
    {
        diagD[i] = 1.0 / (2 * i + 1);
        vectorF[i] = i + 1;

        int baseIndexI = m - i;

        for (int j = max(0, i - m); j < i; ++j)
        {
            int indexIJ = baseIndexI + j;
            matrixAL[i][indexIJ] = 1.0 / (i + j + 1);
        }
    }
}

void SLAUSolverLDLT::returnMatixAfterHilbert()
{
    for (int i = 1; i < n; ++i)
    {
        diagD[i] = 1.0 / (2 * i + 1);
        int baseIndexI = m - i;

        for (int j = max(0, i - m); j < i; ++j)
        {
            int indexIJ = baseIndexI + j;
            matrixAL[i][indexIJ] = 1.0 / (i + j + 1);
        }
    }
}
