#region imports
from NumericalMethods import Secant
from math import cos
#endregion
"""defining functions so we can call them into secant function"""
def fn1(x):
    return x - 3 * cos(x)

def fn2(x):

    return cos(2*x) * (x**3)

"""as stated below, fn1(x)=x-3cos(x) & fn2(x)=cos(2x)x^3 will be implemented in the
secant function in order to find and print roots of fcn(x) (secant function)"""
def main():
    """
    We want to demonstrate the Secant function with:
      1) fn1(x) = x - 3cos(x), x0=1, x1=2, maxiter=5,  xtol=1e-4
      2) fn2(x) = cos(2x)*x^3, x0=1, x1=2, maxiter=15, xtol=1e-8
      3) fn2(x) again with x0=1, x1=2, maxiter=3,  xtol=1e-8
    """
    # 1) fn1
    r1 = Secant(fn1, 1, 2, maxiter=5, xtol=1e-4) #variables set for secant function with x & iter inputs which find root
    print("Root of fn1(x)=x-3cos(x) => {0:.6f}, found in {1} iterations".format(r1[0], r1[1]))

    # 2) fn2 (maxiter=15, xtol=1e-8)
    r2 = Secant(fn2, 1, 2, maxiter=15, xtol=1e-8)
    print("Root of fn2(x)=cos(2x)*x^3 => {0:.6f}, found in {1} iterations (maxiter=15)".format(r2[0], r2[1]))

    # 3) fn2 again (maxiter=3, xtol=1e-8)
    r3 = Secant(fn2, 1, 2, maxiter=3, xtol=1e-8)
    print("Root of fn2(x)=cos(2x)*x^3 => {0:.6f}, found in {1} iterations (maxiter=3)".format(r3[0], r3[1]))

if __name__ == "__main__":
    main()
