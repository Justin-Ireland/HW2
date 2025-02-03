# region imports
import copy as cp
import math
# endregion

# region functions
def GetMinorMatrix(A, row, col):
    '''
    This removes a row and col from a deep copy of A
    :param A: original matrix
    :param row: row index to remove
    :param col: column index to remove
    :return: modified matrix
    '''
    B=cp.deepcopy(A) #make a deep copy of A
    B.pop(row) #removes row from A
    for r in B:  #removes column col from A
        r.pop(col)
    return B

def Determinant(A):
    '''
    Compute the determinant af a nXn matrix using 1st row for cofactors.
    Note: Determinant calls itself if the matrix is larger than 1x1.
    :param A:
    :return:
    '''
    det=0
    if len(A) == 1: #smallest possible minor matrix
        return A[0][0]
    else:
        for c in range(len(A[0])):  #I'm using the first row for cofactors
            sign=(-1)**c  #calculate the sign as prescribed by determinant using cofactors and minors method
            if(A[0][c]!=0): #if cofactor is zero, no contribution to determinant
                det += sign*A[0][c]*Determinant(GetMinorMatrix(A,0,c)) #RECURSIVELY solve determinant
    return det

def main():
    '''
    This main function is to teach list comprehensions and build matrices and test the functions that operate on matrices
    '''
    # First a list comprehension to build a simple list where the element value is the square of its index.
    a = [round(math.sin(x/6*math.pi),4) for x in range(10)]
    # Now, let's build a matrix
    b = [[round(math.sin(x/6*math.pi)*(y+1),4) for x in range(10)] for y in range(10)]
    # Determinant of a matrix
    A=[[1,-1,2],
       [2,3,-3],
       [4,5,1]] #an easy matrix to verify Determinant function. Can you do it by hand?  Can you do it in your head?

    detA=Determinant(A)
    print("A=")
    for r in A:
        print(r)
    print("det(A)=",detA)
    print("")

    AnotherA = [[1, -10, 2, 4],
                [3, 1, 4, 12],
                [9, 2, 3, 4],
                [-1, 2, 7, 3]]

    detAnotherA=Determinant(AnotherA)
    print("A=")
    for r in AnotherA:
        print(r)
    print("det(AnotherA)=",detAnotherA)
    print("")
# endregion

if __name__=='__main__':
    main()
