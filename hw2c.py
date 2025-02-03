from copy import deepcopy
from NumericalMethods import GaussSeidel

def main():
    """This function uses the Gauss-Seidel method to estimate the solution
    to set of N linear equations expressed in matrix form as Ax = b. sets variables for two matrices, solves
    matrices and print output to user"""

    A1 = [[3.0,  1.0, -1.0,  2.0], #matrix 1 3x3
        [1.0,  4.0,  1.0, 12.0],
        [2.0,  1.0,  2.0, 10.0]]

    x_init = [0.0, 0.0, 0.0]  ##starting point set 0 in [x,y,z] form to right of matrix
    solution = GaussSeidel(A1, x_init, Niter=15)
    print("Solution to 3x3 system [3,1,-1; 1,4,1; 2,1,2]:", solution)

    """Now we solve the 4x4 matrix using Gauss Seidel with a starting point [x,y,z,w] ==> zero"""
    A2 = [[ 1.0, -10.0,  2.0,  4.0,  2.0], #matrix 2 4x4
        [ 3.0,   1.0,  4.0, 12.0, 12.0],
        [ 9.0,   2.0,  3.0,  4.0, 21.0],
        [-1.0,   2.0,  7.0,  3.0, 37.0]]


    x_init2 = [0.0, 0.0, 0.0, 0.0]  # initial starting point for [x,y,z,w]
    solution2 = GaussSeidel(A2, x_init2, Niter=15) #call G-S for A2 w/ 15 iterations
    print("Solution to 4x4 system [1,-10,2,4; 3,1,4,12; 9,2,3,4; -1,2,7,3]:", solution2)
    # calls G-S and uses A2 matrix and starting point to solve

    """(blanket statement: Cite ChatGPT to correct code into properly running program for NM.py, hw2a,b,c.py"""
    #endregion

if __name__ == "__main__":
    main()
