#region imports
import Gauss_Elim as GE
from math import sqrt, pi, exp
#endregion

def GPDF(args):
    """Here is where I will define the Gaussian probability density function.
    This requires knowing the population mean and standard deviation.
    To compute the GPDF at any value of x, I just need to compute as stated
    in the homework assignment.
    Step 1:  unpack the args tuple into variables called: x, mu, stDev
    Step 2:  compute GPDF value at x
    Step 3:  return value
    :param args: (x, mean, standard deviation)  tuple in that order
    :return: value of GPDF at the desired x

    The formula is:
        (1 / (sig * sqrt(2*pi))) * exp(-0.5 * ((x - mu)/sig)^2)
    :param args: (x, mu, sig)
    :return: float value of the PDF at x
    """
    x, mu, sig = args
    return (1.0 / (sig * sqrt(2.0 * pi))) * exp(-0.5 * ((x - mu) / sig)**2)


def Simpson(fn, args, N=100):
    """
        This executes the Simpson 1/3 rule for numerical integration
        1. divide the range from x=lhl to x=rhl into an even number of parts.
        2. compute fn at each x value between lhl and rhl
        3. sum the even and odd values of fx as prescribed
        4. return the area beneath the function fx
        :param fn: some function of x to integrate
        :param args: a tuple containing (mean, stDev, lhl, rhl)
        :return: the area beneath the function between lhl and rhl as the approximate integration

        Simpson's from txtbk:h/3[f(a)+f(b)+4(summ.f(x))odd +2(summ.f(x))even]"""
    mu, sig, a, b = args #similar to probability but with lhl=a & rhl=b

    # #check if N is even, if not, increase by 1
    if N % 2 != 0:
        N += 1

    h = (b - a) / N #(rhl-lhl)/100 for area of simpson's
    total = fn((a, mu, sig)) #set var to implement a,mu,sig to fn

    # sum over intermediate points
    for i in range(1, N): #used to find x for each interval
        x_i = a + i * h
        if i % 2 == 0: #if even increment interval and sum with 2fn(x,mu,sig)
            total += 2.0 * fn((x_i, mu, sig))
        else:
            total += 4.0 * fn((x_i, mu, sig)) #if odd increment & sum w/ 4fn(x,mu,sig)

    total += fn((b, mu, sig)) #end rh limit incremented to N=100

    return (h / 3.0) * total


def Probability(PDF, args, c, GT=True):
    """ This is the function to calculate the probability that x is >c or <c depending
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
    :return: probability value"""
    mu, sig = args #unpacked args
    lhl = mu - 5.0 * sig #left hand limit set
    rhl = mu + 5.0 * sig #right hand limit set
    args1 = args, lhl, rhl  # wrap above variables up as args1
    p = Simpson(PDF, (mu, sig, a, b))  # call simpson with PDF and args 1 plugin
    if GT: #if GT=True prob x>c
        return 1.0 - p
    else: #GT=False prob x<c
        return p


def Secant(fcn, x0, x1, maxiter=10, xtol=1e-5):
    """
    Uses the Secant method to find a root of equation fcn(x) = 0 near x0, x1.
    written equation in a form
    fcn = 0 such that when the correct value of x is selected, the fcn actually equals zero (or very close to it).
    :param fcn: function for which we want a root
    :param x0: initial guess
    :param x1: second guess
    :param maxiter: maximum iterations
    :param xtol: tolerance on consecutive x-values
    :return: tuple with: (the final estimate of the root (most recent value of x), number of iterations)
    """
    f0 = fcn(x0)
    f1 = fcn(x1)

    for i in range(maxiter): #loop until maxiter=10
        denom = (f1 - f0) #assign difference to denom
        if abs(denom) < 1e-15: #almost = zero
            # can't proceed if denominator is ~0
            return (x1, i) #stops if denominator is too small
        x2 = x1 - f1 * (x1 - x0) / denom ##secant formula assigned to x2

        if abs(x2 - x1) < xtol: #becomes x1,x2 so function is able to iterate and converge
            return (x2, i + 1)

        """the secant method uses iterated approximations, so we need to
        shift Xs in loop to continue approximations x1,x2 as most recent"""
        ##becomes x1,x2 so function is able to iterate and converge
        x0, x1 = x1, x2
        f0, f1 = f1, fcn(x2)

    return (x2, maxiter)


def GaussSeidel(Aaug, x, Niter=15):
    """
    Solves A x = b via Gauss-Seidel iteration on the augmented matrix Aaug = [A|b].
    Ensures diagonal dominance by reordering rows if possible.  :param Aaug: The augmented matrix from Ax=b -> [A|b]
    :param x:  An initial guess for the x vector. if A is nxn, x is nx1
    :param Niter:  Number of iterations to run the GS method
    :return: the solution vector x
    """
    #reorders rows to make them diag dominant
    Aaug = GE.MakeDiagDom(Aaug)

    n = len(Aaug) #length matrix
    for _ in range(Niter): #loop 15 times
        for i in range(n):
            # sum over j if not equal to i (!=)
            s = 0.0 #starting point
            for j in range(n):
                if j != i: #if j not equal to i, iterate s
                    s += Aaug[i][j] * x[j]
            x[i] = (Aaug[i][n] - s) / Aaug[i][i] #solution vector after iterations

    return x


def main():
    """asked chatgpt to provide blocks to be able to test each numerical method. Here we can test results of GPDF,
    Secant, Gauss Seidel, and Probability"""
    # Optional internal test
    print("NumericalMethods main() quick tests...")

    #GPDF test
    val = GPDF((0, 0, 1))
    print("GPDF(0|mu=0, sig=1) =", val)

    # Probability test: P(x<0) or P(x>0) for N(0,1)
    p_less = Probability(GPDF, (0,1), 0, GT=False)
    p_greater = Probability(GPDF, (0,1), 0, GT=True)
    print("P(x<0|N(0,1)) =", p_less)
    print("P(x>0|N(0,1)) =", p_greater)

    # Quick Secant example: f(x)=x^2-2
    def f(x): return x*x - 2
    root, iters = Secant(f, 1, 2)
    print("Root of x^2-2 ~", root, "in", iters, "iterations")

    # Quick GaussSeidel example
    # 2x+y=5, x+3y=9 in matrix form [[2,1,5],[1,3,9]]
    Aaug_test = [[2.0,1.0,5.0],[1.0,3.0,9.0]]
    guess = [0.0,0.0]
    sol = GaussSeidel(Aaug_test, guess, Niter=15)
    print("GaussSeidel solution for 2x+y=5, x+3y=9 =>", sol)


if __name__ == "__main__":
    main()
