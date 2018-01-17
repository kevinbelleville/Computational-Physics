##
## Consider the quantum anharmonic oscillator defined by:
##              H_lambda = H_0 + lambda*x^4
##
## Step one: Determine the matrix elements of H_lambda in the harmonic
## oscillator basis, B_0 using raising and lowering operators.
##      raising operator: A_raise
##      lowering operator: A_lower
##
## Separate the two quantities, H_0 and lambda*x^4
##
## H_0 = A_raise*A_lower + 1/2*I
## lambda*x^4 = lambda*(A_raise + A_lower)^4
##
## H_lambda = H_0 + lambda*x^4
## H_lambda = A_raise*A_lower + 1/2*I + lambda*(A_raise + A_lower)^4
##
## To find the matrix elements, we apply bra <m| and ket |n> to the
## quantities with the raising the lowering operators. The delta function
## is the dirac delta function.
##
## <m| H_0 |n> = (n + 1/2) delta(m,n)
## <m| lambda*x^4 |n> = lambda * <m| x^4 |n>
## lambda * <m| x^4 |n> = lambda * [ (1/4)*(6n^2 + 6n + 3) delta(m,n)
##                          (1/2)*sqrt((n+1)(n+2)) delta]

import math
import numpy as np
import matplotlib.pyplot as plt

def A_raise(n):
    """ Raises an operator. Returns the value, and the raised n. """
    value = math.sqrt(n+1)
    new_n = n+1
    return value, new_n

def A_lower(n):
    """ Lowers an operator. Returns the value, and the lowered n. """
    value = math.sqrt(n)
    new_n = n-1
    return value, new_n

def delta(m, n):
    """ Computes the dirac delta function for m and n. If m == n, return 1,
        otherwise return 0. """

## Let's test this by calulating  <m| H_0 |n>
## H_0 = A_raise * A_lower + (1/2)I

def calc_H_0(m, n):
    first = A_lower(n)
    constants = first[0]
    second = A_raise(first[1])
    constants *= second[0]
    delta()
