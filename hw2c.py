from copy import deepcopy

def main():
    pass

if __name__=="__main":
    main()


# Reed's portion

from NumericalMethods import GaussSeidel


# Implement Gauss-Seidel Method without numpy and without b/b2 vectors
def GaussSeidel(A, x_init, max_iter=100, tol=1e-10):
    n = len(A)
    x = x_init.copy()

    for _ in range(max_iter):
        x_new = x.copy()
        for i in range(n):
            sigma = sum(A[i][j] * x_new[j] for j in range(n) if j != i)
            x_new[i] = (A[i][n] - sigma) / A[i][i]

        # Check for convergence (using list-based norm calculation)
        norm = sum((x_new[i] - x[i]) ** 2 for i in range(n)) ** 0.5
        if norm < tol:
            break
        x = x_new

    return x


def main():
    # Augmented coefficient matrix (A with last column as the constants)
    A = [[3, 1, -1, 2],
        [1, 4, 1, 12],
        [2, 1, 2, 10]]

    # Initial guess (optional)
    x_init = [0, 0, 0]  # Starting guess for [x, y, z]

    # Solve the system using Gauss-Seidel method
    solution = GaussSeidel(A, x_init)

    print("Solution estimated by Gauss-Seidel:", solution)

    # Another example for a different system
    A2 = [[1, -10, 2, 4, 2],
        [3, 1, 4, 12, 12],
        [9, 2, 3, 4, 21],
        [-1, 2, 7, 3, 37]]

    # Initial guess (optional)
    x_init2 = [0, 0, 0, 0]  # Starting guess for [x, y, z, w]

    # Solve the second system using Gauss-Seidel method
    solution2 = GaussSeidel(A2, x_init2)

    print("Solution estimated by Gauss-Seidel for the second system:", solution2)


if __name__ == "__main__":
    main()

# End Reed's portion