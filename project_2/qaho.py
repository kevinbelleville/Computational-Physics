import numpy as np
from numpy.linalg import eig, eigh
import matplotlib.pyplot as plt
from numpy.polynomial.hermite import hermval

np.set_printoptions(precision=3, suppress=True)

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

## variables
N = 100
lam = 1

## generate matrices with zeros
H_0 = np.zeros([N+1, N+1])
H = np.zeros([N+1, N+1])
H_lam = np.zeros([N+1, N+1])

## populate the matrices with the right values
## for H_0 = (n + 0.5) delta(m,n)
for m in range(N+1):
    for n in range(N+1):
        if m == n:
            H_0[m,n] = n + 0.5

## for H = (n + 0.5) delta(m, n) + lambda*sqrt(n+1)*delta(m,n+1)
##                              + lambda*sqrt(n)*delta(m,n-1)
for m in range(N+1):
    for n in range(N+1):
        if m == n:
            H[m, n] = n + 0.5
        elif m == n + 1:
            H[m, n] = lam * np.sqrt(n + 1)
        elif m == n - 1:
            H[m, n] = lam * np.sqrt(n)

## populate H_lam...
for m in range(N+1):
    for n in range(N+1):
        if m == n:
            H_lam[m, n] = n + 0.5 + lam*(0.25)*( 6 * n**2 + 6 * n + 3)
        elif (m == (n+2)):
            H_lam[m, n] = lam*(n + (3/2)) * np.sqrt((n+1)*(n+2))
        elif (m == (n-2)):
            H_lam[m, n] = lam*(n - (1/2)) * np.sqrt((n-1)*n)
        elif (m == (n+4)):
            H_lam[m, n] = lam*(1/4) * np.sqrt((n+1)*(n+2)*(n+3)*(n+4))
        elif (m == (n-4)):
            H_lam[m, n] = lam*(1/4) * np.sqrt((n-3) * (n-2) * (n-1) * n)


## just the matrix reps
# print(H_0)
# print(eig(H_0)[0])
# print(H)
# print(H_lam)
# print(eigh(H_lam)[0])


def generate_H_lam_evals(lam, N):
    H_lam = np.zeros([N+1, N+1])
    for m in range(N+1):
        for n in range(N+1):
            if m == n:
                H_lam[m, n] = n + 0.5 + lam*(0.25)*( 6 * n**2 + 6 * n + 3)
            elif (m == (n+2)):
                H_lam[m, n] = lam*(n + (3/2)) * np.sqrt((n+1)*(n+2))
            elif (m == (n-2)):
                H_lam[m, n] = lam*(n - (1/2)) * np.sqrt((n-1)*n)
            elif (m == (n+4)):
                H_lam[m, n] = lam*(1/4) * np.sqrt((n+1)*(n+2)*(n+3)*(n+4))
            elif (m == (n-4)):
                H_lam[m, n] = lam*(1/4) * np.sqrt((n-3) * (n-2) * (n-1) * n)
    return eigh(H_lam)[0]

#print(generate_H_lam_evals(0, 100)[0:4])
#print(generate_H_lam_evals(0.2, 100)[0:4])
#print(generate_H_lam_evals(0.4, 100)[0:4])
#print(generate_H_lam_evals(0.6, 100)[0:4])
#print(generate_H_lam_evals(0.8, 100)[0:4])
#print(generate_H_lam_evals(1, 100)[0:4])

##
## Plot first four energy levels over the range(0 <= lam <= 1)
##

lam_vals = []
E_vals_0 = []
E_vals_1 = []
E_vals_2 = []
E_vals_3 = []
E_vals_4 = []
E_vals_5 = []
E_vals_6 = []
basis_size = 100

for i in range(0, 21):
    l = i*(1/20)
    lam_vals.append(l)
    E_vals_0.append(generate_H_lam_evals(l, basis_size)[0])
    E_vals_1.append(generate_H_lam_evals(l, basis_size)[1])
    E_vals_2.append(generate_H_lam_evals(l, basis_size)[2])
    E_vals_3.append(generate_H_lam_evals(l, basis_size)[3])
    E_vals_4.append(generate_H_lam_evals(l, basis_size)[4])
    E_vals_5.append(generate_H_lam_evals(l, basis_size)[5])
    E_vals_6.append(generate_H_lam_evals(l, basis_size)[6])


plt.plot(lam_vals, E_vals_0, "c-")
plt.plot(lam_vals, E_vals_1, "r-")
plt.plot(lam_vals, E_vals_2, "g-")
plt.plot(lam_vals, E_vals_3, "y-")
plt.plot(lam_vals, E_vals_4, "b-")
plt.plot(lam_vals, E_vals_5, "m-")
plt.plot(lam_vals, E_vals_6, "k-")
#plt.show()

## list of e_vals
list_E_vals = [E_vals_0, E_vals_1, E_vals_2,
                E_vals_3, E_vals_4, E_vals_5, E_vals_6]

## differences between E levels
def delta_E(E_vals_i, E_vals_i_1):
    delta = []
    for ind, val in enumerate(E_vals_i):
        delta.append(E_vals_i_1[ind] - val)
    return delta

delta_1_0 = delta_E(E_vals_0, E_vals_1)
plt.plot(lam_vals, delta_1_0, "r-")
#plt.show()

colors = ["c-", "r-", "g-", "y-", "b-", "m-", "k-"]

for ind, evals in enumerate(list_E_vals):
    if ind == len(list_E_vals) - 1:
        break
    plt.plot(lam_vals, delta_E(evals, list_E_vals[ind + 1]), colors[ind])
#plt.show()


## hermval stuff
## hermval accepts the parameters:
## c, the length of n+1
## x, the point at which to calculate
## returns an ndarray
# print(hermval([0,1,2,3,4], [0, 0.2, 0.4, 0.6, 0.8, 1]))


## convergence of the method wrt the basis size N
## we shall use E_0( lam = 1 )
base_size = []
for i in range(5, 101):
    base_size.append(i)

inv_base_size = []
for base in base_size:
    inv_base_size.append(1/base)

base_e = []
for base in base_size:
    base_e.append(generate_H_lam_evals(1, base)[0])

plt.plot(base_size, base_e, "r-", base_size, base_e, "r.")
#plt.show()

#print(base_e)

print(eigh(H_lam)[1][:,0])


## eig function
# print(eig(H_0)[0])
# print(eig(H)[0])

#
# print(eig(H_0)[1][:,0])
# print(eig(H)[1][:,0])

#
# v_0 = eig(H)[1][:,0]
# eval_0 = eig(H)[0][0]
# print(np.dot(H, v_0))
# print(eval_0 * v_0)

##
## Wavefunctions
##
## it is possible to think about the vectors |n> in the basis B_0
## as wavefunctions instead of abstract, and in factor the wavefunction
## psi_n that corresponds to |n> is:
##
## psi_n(x) = [(2^n)*(n!)*(sqrt(pi))]^(-1/2) * e^(-x^2/2) * h_n(x)
##
## h_0(x) = 1
## h_1(x) = 2x
## h_n+1(x) = 2x*h_n(x) - 2n*h_(n-1)(x)
##

## create memoization method

def memoize(f):
    memo = {}
    def helper(*args):
        if args not in memo:
            memo[args] = f(*args)
        return memo[args]
    return helper


## create memoized Hermite polynomial

@memoize
def h(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return 2*x
    else:
        return 2*x*h(n - 1, x) - 2*(n - 1)*h(n-2, x)

## create facotrial method for psi

@memoize
def factorial(x):
    if (x == 0) or (x == 1):
        return 1
    else:
        return x * factorial(x - 1)



## create psi method

def psi(n, x):
    first = (2**n)*(factorial(n))*(np.sqrt(np.pi))
    first = first**(-1/2)
    second = np.e**(-(x**2)/2)
    herm = h(n, x)
    soln = first*second*herm
    return soln

## create new psi method for anharmonic version
def apsi(n, x, l):
    hpsi = psi(n,x)
    correction = (l/4)*(6*(n**2) + 6*n + 3)/(2*n + 1)
    psi_first_order_corr = hpsi*correction
    soln = hpsi + (l * psi_first_order_corr)
    return soln


x_vals = []
psi_0_vals = []
psi_1_vals = []
psi_2_vals = []
psi_3_vals = []
for i in range(101):
    x = -4+(2/25)*i
    x_vals.append(x)
    psi_0_vals.append(apsi(0,x,1))
    psi_1_vals.append(apsi(1,x,1))
    psi_2_vals.append(apsi(2,x,1))
    psi_3_vals.append(apsi(3,x,1))

plt.gcf().clear()
plt.plot(x_vals, psi_0_vals, "r-", label="Zeroth")
plt.plot(x_vals, psi_1_vals, "c-", label="First")
plt.plot(x_vals, psi_2_vals, "g-", label="Second")
plt.plot(x_vals, psi_3_vals, "m-", label="Third")

plt.legend()

plt.axhline(0, color='k', linewidth=1)
plt.axvline(0, color='k', linewidth=1)

plt.show()
