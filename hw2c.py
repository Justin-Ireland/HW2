from copy import deepcopy

"""This function uses the Gauss-Seidel method to estimate the solution 
to set of N linear equations expressed in matrix form as Ax = b.
A convergence check is implemented and uses changes between iterations to check convergence."""

def GaussSeidel(A, x_init, max_iter=100, tol=1e-10):
    n = len(A) #size of matrix
    x = x_init.copy()

    for _ in range(max_iter): #loop to 100 iterations maximum
        x_new = x.copy() #previous iteration
        for i in range(n): #loops through each value in A
            sigma = sum(A[i][j] * x_new[j] for j in range(n) if j != i) #sums if j not equal to i
            x_new[i] = (A[i][n] - sigma) / A[i][i] #gauss-seidel formula

        # Check for convergence (using list-based norm calculation)
        norm = sum((x_new[i] - x[i]) ** 2 for i in range(n)) ** 0.5
        if norm < tol:
            break
        x = x_new

    return x

"""The main function sets variables for two matrices and will use the above Gauss Seidel method to solve
 matrices and print output to user"""

def main():
    # Augmented coefficient matrix (A with last column as the constants)
    A = [[3, 1, -1, 2],
        [1, 4, 1, 12],
        [2, 1, 2, 10]]

    # Initial guess
    x_init = [0, 0, 0]  #starting point set 0 in [x,y,z] form to right of matrix

    # Solve the system using Gauss-Seidel method
    solution = GaussSeidel(A, x_init)

    print("Solution estimated by Gauss-Seidel:", solution)

    # Another example for a different system
    A2 = [[1, -10, 2, 4, 2],
        [3, 1, 4, 12, 12],
        [9, 2, 3, 4, 21],
        [-1, 2, 7, 3, 37]]

    # Initial guess (optional)
    x_init2 = [0, 0, 0, 0]  # Starting point for [x, y, z, w] in 4x5

    # Solve the second system using Gauss-Seidel method
    solution2 = GaussSeidel(A2, x_init2) # calls G-S and uses A2 matrix and starting point to solve

    print("Solution estimated by Gauss-Seidel for the second system:", solution2)

if __name__=="__main":
    main()