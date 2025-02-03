#region imports
from NumericalMethods import Secant
from math import cos
#endregion

#region function definitions
def fn1(x):
    pass

def fn2(x):
    pass

def main():
    """
       fn1:  x-3cos(x)=0; with x0=1, x1= 2, maxiter = 5 and xtol = 1e-4
       fn2:  cos(2x)*x**3; with x0=1, x1= 2, maxiter = 15 and xtol = 1e-8
       fn2:   with x0=1, x1= 2, maxiter = 3 and xtol = 1e-8

       I observe that for functions 2 and 3, the answer should be pi/2 or about 1.57
    :return: nothing, just print results
    """
    r1 = Secant(fn1, 1, 2, 5,1e-4)
    r2 = Secant(fn2, 1,2,15, 1e-8)
    r3 = Secant(fn2,1,2,3,1e-8)
    print("root of fn1 = {root:0.4f}, after {iter :0d} iterations".format(root=r1[0], iter=r1[1]))
    #etc.
    pass
#endregion

if __name__=="__main__":
    main()




#Reed's portion

#region imports
from NumericalMethods import Secant
from math import cos
#endregion

def fn1(x):
    """
    fn1(x) = x - 3*cos(x).
    We'll try to find a root near [1, 2].
    """
    return x - 3*cos(x)

def fn2(x):
    """
    fn2(x) = cos(2x) * x^3
    We'll attempt to find a root near [1, 2].
    """
    return cos(2*x) * (x**3)

def main():
    """
    We want to demonstrate the Secant function with:
      1) fn1(x) = x - 3cos(x), x0=1, x1=2, maxiter=5,  xtol=1e-4
      2) fn2(x) = cos(2x)*x^3, x0=1, x1=2, maxiter=15, xtol=1e-8
      3) fn2(x) again with x0=1, x1=2, maxiter=3,  xtol=1e-8
    """
    # 1) fn1
    r1 = Secant(fn1, 1, 2, maxiter=5, xtol=1e-4)
    print("Root of fn1(x)=x-3cos(x) => {0:.6f}, found in {1} iterations".format(r1[0], r1[1]))

    # 2) fn2 (maxiter=15, xtol=1e-8)
    r2 = Secant(fn2, 1, 2, maxiter=15, xtol=1e-8)
    print("Root of fn2(x)=cos(2x)*x^3 => {0:.6f}, found in {1} iterations (maxiter=15)".format(r2[0], r2[1]))

    # 3) fn2 again (maxiter=3, xtol=1e-8)
    r3 = Secant(fn2, 1, 2, maxiter=3, xtol=1e-8)
    print("Root of fn2(x)=cos(2x)*x^3 => {0:.6f}, found in {1} iterations (maxiter=3)".format(r3[0], r3[1]))

if __name__ == "__main__":
    main()

#End Reed's portion