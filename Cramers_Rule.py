#region imports
from Determinant import Determinant
from Gauss_Elim import popColumn, insertColumn
from copy import deepcopy
#endregion
def Cramer(A, b):
    '''
    See lecture notes for Cramer's method to solve the matrix equation [A][x]=[b]
    A is a matrix
    b is a column vector
    x is the unknow vector for which we are trying to solve

    Step 1:  compute determinant of A
    Step 2:  replace first column of A with b vector, call this matrix AA
    Step 3:  compute first row of x (a single value) by x[i]=det[AA]/det[A]
    Step 4:  repeat steps 3 & 4 for the remaining elements of x
    '''
    n = len(A)
    x = [0]*n
    detA = Determinant(A)
    for colIndex in range(n):  # generate AA for each col
        AA=popColumn(A,colIndex)[1]
        AA=insertColumn(AA,b, colIndex)
        detAA = Determinant(AA)  # get the determinant of AA
        x[colIndex] = detAA / detA  # the essence of Cramer's Rule
    return x

def main():
    '''
    This is used to explore Cramer's method
    '''
    A = [[1, -2, 3, 4], [5, 6, 7, 8], [-9, 10, -11, 6], [5, 4, -3, 2]]
    b = [1, 2, 3, 4]

    print(f"detA={Determinant(A):0.3f}")
    print("x=", Cramer(A, b))

    A = [[-5, 1, -5, 0, 1, -4], [5, 0, 3, 5, 3, 5], [-2, -2, 1, 4, 3, -5],
                  [4, 5, 0, 3, 4, -1], [-5, -2, -5, 5, -2, -2], [4, 5, 5, 0, 0, -2]]
    b = [-95, -45, 49, -50, 90, 30]

    print(f"detA={Determinant(A):0.3f}")
    print("x=", Cramer(A, b))

if __name__=='__main__':
    main()
