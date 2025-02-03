#region imports
from math import sqrt, pi, exp
from NumericalMethods import GPDF, args, Simpson, Probability

#endregion

#region function definitions
def main():
    """
    I want to integrate the Gaussian probability density function between
    a left hand limit = (mean - 5*stDev) to a right hand limit = (c).  Here
    is my step-by-step plan:
    1. Decide mean, stDev, and c and if I want P(x>c) or P(x<c).
    2. Define args tuple and c to be passed to Probability from NM
    3. Pass args, and a callback function (GPDF) to Probability
    4. In probability, pass along GPDF to Simpson along with the appropriate args tuple
    5. Return the required probability from Probability and print to screen.
    :return: print results to screen.
    """
    #region testing user input
    # The following code solicits user input through the CLI.
    #float changes user input to number
    mean = float(input("Population mean? "))
    stDev = float(input("Standard deviation?"))
    c = float(input("c value?"))
    GT = True if input("Probability greater than c?").lower() in ["y","yes","true"] else "False"
    #no need to float() GT b/c boolean

    """following code calls probability function and applies the user inputs above
    to test if x is above or below rhl c and print probability based of range x"""

    probx = Probability(GPDF, args, c,GT) #assign variable to probability function to condense it
    if GT: #if GT true print P(x>c | N(m,stdev)):probability
        print(f"P(x>{c}|N({mean},{stDev}))={probx:.4f}")
    else: #GT False then print P(x<c | N(m,stdev)):probability
        print(f"P(x<{c}|N({mean},{stDev}))={probx:.4f}")
    #endregion
#endregion

if __name__ == "__main__":
    main()