# NSTU 5th Semester Numerical Methods Practical Task 1

# LDL<sup>T</sup> Decomposition Project

This project provides an implementation of the LDL<sup>T</sup> decomposition for symmetric matrices. It includes utilities for matrix operations and solving linear systems based on the LDL<sup>T</sup> decomposition method. The decomposition splits a matrix \( A \) into three components: \( L \) (lower triangular matrix), \( D \) (diagonal matrix with entries ±1), and \( L^T \) (transpose of \( L \)).

## Solving the Linear System \( Ax = F \) After LDL<sup>T</sup> Decomposition

After performing LDL<sup>T</sup> decomposition, where the matrix \( A \) is decomposed into \( L \) (lower triangular matrix), \( D \) (diagonal matrix), and \( L^T \) (transpose of \( L \)), solving the linear system \( Ax = F \) is carried out in the following steps:

### Steps to Solve \( Ax = F \)

1. **Decomposition**: The matrix \( A \) is decomposed into \( LDL^T \):
   \[
   A = L D L^T
   \]
   where:
   - \( L \) is a lower triangular matrix.
   - \( D \) is a diagonal matrix with entries ±1.
   - \( L^T \) is the transpose of \( L \).

2. **Transform the System**: Substitute \( A = LDL^T \) into the linear system \( Ax = F \):
   \[
   LDL^T x = F
   \]
   Let \( y = L^T x \). This transforms the system into:
   \[
   LD y = F
   \]

3. **Solve for \( y \)**:
   - **Forward Substitution**: Solve the lower triangular system \( LD y = F \) for \( y \).
     \[
     L y = z
     \]
     where \( z \) is an intermediate vector obtained from:
     \[
     z = F
     \]
     Forward substitution is used to solve \( L y = z \).

4. **Solve for \( x \)**:
   - **Diagonal Scaling**: Solve the diagonal system \( D y = z \):
     \[
     y_i = \frac{z_i}{D_i}
     \]
     where \( D_i \) is the diagonal element.
   - **Back Substitution**: Solve the upper triangular system \( L^T x = y \):
     \[
     L^T x = y
     \]
     Back substitution is used to solve \( L^T x = y \).

## Project Structure

### Data Directory

- **`AL.txt`**: File containing the lower triangular part of the matrix \( AL \).
- **`D.txt`**: File containing the diagonal matrix \( D \).
- **`F.txt`**: File containing the vector \( F \) for solving the system \( Ax = F \).
- **`input.txt`**: File containing the matrix dimensions \( n \) and \( k \).
- **`x.txt`**: Output file for the solution vector \( x \).

## How to Build and Run

1. **Build the Project**: Use the `Makefile` to compile the project. The build targets will create executables for different floating point precisions.
   ```sh
   make
