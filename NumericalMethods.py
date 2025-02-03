#region imports
from openpyxl.styles.builtins import total

import Gauss_Elim as GE  # this is the module from lecture 2 that has useful matrix manipulation functions
from math import sqrt, pi, exp
#endregion

#region function definitions
def Probability(PDF, args, c, GT=True):
    """
    This is the function to calculate the probability that x is >c or <c depending
    on the GT boolean.
    Step 1:  unpack args into mu and stDev
    Step 2:  compute lhl and rhl for Simpson
    Step 3:  package new tuple args1=(mu, stDev, lhl, rhl) to be passed to Simpson
    Step 4:  call Simpson with GNPDF and args1
    Step 5:  return probability
    :param GPDF: the probability density function to be integrated
    :param args: a tuple with (mean, standard deviation)
    :param c: value for which we ask the probability question
    :param GT: boolean deciding if we want probability x>c (True) or x<c (False)
    :return: probability value
    """
    mu, sig = args #unpacked args
    lhl= mu - 5 * sig #left limit set
    rhl= c = mu + 5 * sig #right limit set, =c
    args1 = args, lhl, rhl #wrap above variables up as args1
    p = Simpson(PDF, args1) #call simpson with PDF and args 1 plugins
    if GT: #if GT=True prob x>c
        return 1.0 - p
    else: #GT=False prob x<c
        return p

def GPDF(args):
    """
    Here is where I will define the Gaussian probability density function.
    This requires knowing the population mean and standard deviation.
    To compute the GPDF at any value of x, I just need to compute as stated
    in the homework assignment.
    Step 1:  unpack the args tuple into variables called: x, mu, stDev
    Step 2:  compute GPDF value at x
    Step 3:  return value
    :param args: (x, mean, standard deviation)  tuple in that order
    :return: value of GPDF at the desired x
    """
    # Step 1: unpack args
    x, mu, sig = args
    # step 2: compute GPDF at x
    fx = (1 / (sig * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sig) ** 2)
    # step 3: return value
    return fx

def Simpson(fn, args, N=100):
    """
    This executes the Simpson 1/3 rule for numerical integration (see page 832, Table 19.4).
    As I recall:
    1. divide the range from x=lhl to x=rhl into an even number of parts. Perhaps 20?
    2. compute fx at each x value between lhl and rhl
    3. sum the even and odd values of fx as prescribed
    4. return the area beneath the function fx
    :param fx: some function of x to integrate
    :param args: a tuple containing (mean, stDev, lhl, rhl)
    :return: the area beneath the function between lhl and rhl as the approximate integration
    """
    """Simpson's from txtbk:h/3[f(a)+f(b)+4(summ.f(x))odd +2(summ.f(x))even]"""
    mu, sig, a, b = args #similar to probability but with lhl=a & rhl=b
    h = (b-a)/N #(rhl-lhl)/100 for area of simpson's
    total = fn((a,mu,sig)) #set var to implement a,mu,sig to fn

    if N % 2 != 0: #check if N is even, if not, increase by 1
        N +=1
    for i in range(1,N): #used to find x for each interval
        x = a + i * h
        if i % 2 ==0: #if even increment interval and sum with 2fn(x,mu,sig)
            total += 2*fn((x,mu,sig))
        else: total += 4*fn((x,mu,sig)) #if odd increment & sum w/ 4fn(x,mu,sig)
    total += fn((b,mu,sig)) #end rh limit incremented to N=100
    return (h/3)*total



def Secant(fcn, x0, x1, maxiter=10, xtol=1e-5):
    """
    This function implements th Secant method to find the root of an equation.  You should write your equation in a form
    fcn = 0 such that when the correct value of x is selected, the fcn actually equals zero (or very close to it).
    :param fcn: the function for which we want to find the root
    :param x0: x value in neighborhood of root (or guess 1)
    :param x1: another x value in neighborhood of root (or guess x0+1)
    :param maxiter: exit if the number of iterations (new x values) equals this number
    :param xtol:  exit if the |xnewest - xprevious| < xtol
    :return: tuple with: (the final estimate of the root (most recent value of x), number of iterations)
    """
    f0 = fcn(x0)
    f1 = fcn(x1)
    for i in range(maxiter): #loop until maxiter=10
        denom =  (f1-f0) #assign difference to denom
        if abs(denom) <1e-15:  #~0
            return (x1, i) #stops if denominator is too small

        x2= x1 - fcn(x1)*(x1-x0)/denom #secant formula assigned to x2

        if abs(x2-x1) < xtol: #check if diff <xtol, if so, converges
            return (x2, i + 1)

        """the secant method uses iterated approximations, so we need to
        shift Xs in loop to continue approximations x1,x2 as most recent"""
        x0, x1 = x1, x2 #becomes x1,x2 so function is able to iterate and converge
        f0, f1 = f1, fcn(x2)



    return (x2, maxiter)


def GaussSeidel(Aaug, x, Niter = 15):
    """
    This should implement the Gauss-Seidel method (see page 860, Tabl 20.2) for solving a system of equations.
    :param Aaug: The augmented matrix from Ax=b -> [A|b]
    :param x:  An initial guess for the x vector. if A is nxn, x is nx1
    :param Niter:  Number of iterations to run the GS method
    :return: the solution vector x
    """
    Aaug = GE.MakeDiagDom(Aaug) #reorders rows to make them diag dominant

    n = len(Aaug) #length matrix
    for _ in range(Niter): #loop 15 times
        for i in range(n):
            s = 0 #starting point
            for j in range(n):
                if j != 1: #if j not equal to one, iterate s
                    s += Aaug[i][j] * x[j] #iter
            x[i] = (Aaug[i][n]-s)/Aaug[i][i] #solution vector after iters
    return x

def main():
    '''
    This is a function I created for testing each numerical method locally.
    :return: None
    '''
    #region testing GPDF
    fx = GPDF((0,0,1)) #GPDF (x,mu,sig)
    print("{:0.5f}".format(fx))  # Does this match the expected value?
    #edregion

    #region testing Simpson
    p=Simpson(GPDF,(0,1,-5,0)) # should return 0.5
    print("p={:0.5f}".format(p))  # Does this match the expected value?
    #endregion

    #region testing Probability
    p1 = Probability(GPDF, (0,1),0,True) #prob if true
    print("p1={:0.5f}".format(p1))  # Does this match the expected value?
    #endregion
    p2 = Probability(GPDF, (0,1),0,False) #prob if false
    print("p2={:0.5f}".format(p2))

#endregion

#region function calls
if __name__ == '__main__':
    main()
#endregion