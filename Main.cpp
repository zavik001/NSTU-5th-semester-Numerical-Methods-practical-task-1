#include "LDLT.hpp"

int main()
{
    try
    {
        string inputFilePath = "data/input.txt";
        string alFilePath = "data/AL.txt";
        string dFilePath = "data/D.txt";
        string fFilePath = "data/F.txt";
        string xFilePath = "data/X.txt";

        LDLT ldlt(inputFilePath, alFilePath, dFilePath, fFilePath, xFilePath);
        ldlt.printRestoredMatrix();
        ldlt.performLDLtDecomposition();
        ldlt.printRestoredMatrix();
        ldlt.solveLinearSystem();
        ldlt.returnMatix();
        ldlt.printRestoredMatrix();
        ldlt.printVectors();
        ldlt.printMultiplyMatrixToVectorX();
    }
    catch (const exception &e)
    {
        cerr << "EROR: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
