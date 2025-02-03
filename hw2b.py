#region imports
from NumericalMethods import Secant
from math import cos
#endregion

#region function definitions
def fn1(x):
    return x - 3*cos(x)

"""as stated below, fn1(x)=x-3cos(x) & fn2(x)=cos(2x)x^3 will be implemented in the
secant function in order to find root of fcn(x)"""

def fn2(x):
    return cos(2*x) * (x**3)

def main():
    """
       fn1:  x-3cos(x)=0; with x0=1, x1= 2, maxiter = 5 and xtol = 1e-4
       fn2:  cos(2x)*x**3; with x0=1, x1= 2, maxiter = 15 and xtol = 1e-8
       fn2:   with x0=1, x1= 2, maxiter = 3 and xtol = 1e-8

       I observe that for functions 2 and 3, the answer should be pi/2 or about 1.57
    just print results
    """
    r1 = Secant(fn1, 1, 2, 5,1e-4)
    r2 = Secant(fn2, 1,2,15, 1e-8)
    r3 = Secant(fn2,1,2,3,1e-8)
    #variables set for secant function with x & iter inputs which find root
    #following code prints results
    print("root of fn1 = {root:0.4f}, after {1} iterations".format(root=r1[0], iter=r1[1]))

    print("root of fn2 = {root:0.4f}, after {1} iterations (maxiter=15)".format(root=r2[0], iter=r2[1]))

    print("root of fn2 = {root:0.4f}, after {1} iterations (maxiter=3)".format(root=r3[0], iter=r3[1]))

#endregion

if __name__=="__main__":
    main()