# Numerical Analysis Assignment 2 code
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.optimize import fminbound
import pandas as pd

# -------------------------------- Question 1: --------------------------------
# g(x) = tan(sin(x^3)) on the interval [0, 2]

# 1. // Compute the interpolating polynomials of degree n = 5 and n = 10 using equally
#    spaced interpolation points //

# SETTING UP THE QUESTION
# function to approximate
def g(x):
    return np.tan(np.sin(x**3))


# FUNCTIONS TO USE IN COMPUTING THE INTERPOLATING POLYNOMIALS - getting matrix equation Ax = b to solve for coeff. x
# Function to get the A matrix
# INPUTS:
#   - n = number of subintervals to split the interval into (such that there are n + 1 points)
#   - p = list of the n + 1 points
# OUTPUT:
#   - The A matrix in Ax = b (to solve for polynomial coefficients x)
def A_matrix(n, points):
    power = np.linspace(n, 0, n + 1)
    return np.array([[p**j for j in power] for p in points])


# Function to get the b vector
# INPUTS:
#   - p = list of the n + 1 points (same p as used in A_matrix(n, p) function)
#   - function = the function that is being approximated by the polynomial interpolation
# OUTPUTS:
#   - the b vector in Ax = b (to solve for polynomial coefficients x)
def b_matrix(points, function):
    return np.array([[function(p)] for p in points])


# Function for creating the nth degree polynomial
# INPUTS:
#   - n = degree of the polynomial
#   - coeff = coefficients starting from highest power of x down all the way to constant
# OUTPUT:
#   - polynomial of degree n with coefficients given, as a function of x
def poly(n, coeff, x):
    power = np.linspace(n, 0, n + 1)
    parts = [coeff[i] * x**power[i] for i in range(n + 1)]
    return sum(parts)



# COMPUTING THE INTERPOLATING POLYNOMIAL FOR n = 5
# n = 5 degree polynomial (6 unknowns) -> n + 1 = 6 points, n intervals
n = 5
nr_points = n + 1
points = np.linspace(0, 2, nr_points)

# A matrix
A = A_matrix(n, points)

# b matrix
b = b_matrix(points, g)

# polynomial coefficients
coeff_n5 = np.linalg.solve(A, b)

# creating dataframe
coeff_n5_table = pd.DataFrame(coeff_n5)
coeff_n5_table.columns = ["coefficients"]
coeff_n5_table.index = ['a5', 'a4', 'a3', 'a2', 'a1', 'a0']

# function for polynomial
def poly_n5(x):
    return poly(5, coeff_n5, x)


# REPEATING FOR THE INTERPOLATING POLYNOMIAL OF DEGREE n = 10
n = 10
nr_points = n + 1
points = np.linspace(0, 2, nr_points)

# A matrix
A = A_matrix(n, points)

# b matrix
b = b_matrix(points, g)

# coefficients
coeff_n10 = np.linalg.solve(A, b)

coeff_n10_table = pd.DataFrame(coeff_n10)
coeff_n10_table.columns = ["coefficients"]
coeff_n10_table.index = ['b10', 'b9', 'b8', 'b7', 'b6', 'b5', 'b4', 'b3', 'b2', 'b1', 'b0']

# n = 10 degree polynomial estimate
def poly_n10(x):
    return poly(10, coeff_n10, x)


# // Plot the interpolants on the same graph //
x = np.linspace(0, 2, 100)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.spines['bottom'].set_position(('data', 0))
ax.spines['left'].set_position(('data', 0))
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.plot(x, poly_n5(x), "b", label="n = 5")
plt.plot(x, poly_n10(x), "r", label="n = 10")
plt.plot(x, g(x))

plt.legend(loc="lower center")
#plt.show()

# 2. // Compute interpolating polynomials of degree n = 5 and n = 10 based on zeros of
#       the appropriate Chebyshev polynomial.
# Procedure: Find zeroes of chebyshev polynomial with degree (n + 1) -> will have n + 1 zeros
#            to be used as the  n + 1 interpolation points for the polynomial approximation of degree n


# FOR THE n = 5 INTERPOLATING POLYNOMIAL WITH POINTS EQUAL TO ZEROS OF m = n + 1 = 6 th degree ChebPoly:
# 5 + 1 = 6th degree Chebyshev polynomial
m = 6

# Calculating zeros of the 6th degree chebyshev polynomial in [-1, 1]
# Function for calculating zeros of mth degree chebyshev polynomial
def cheb_zeros(m):
    zeros = []
    for i in range(m):
        zeros.append(np.cos((2*i + 1) * np.pi/(2*m)))
    return zeros

zeros = cheb_zeros(m)

# transforming zeros to interval [a, b] = [0, 2]
def tfmd_zeros(a, b, zeros):
    return [0.5 * (a + b) + 0.5 * (b - a)*p for p in zeros]

zeros = tfmd_zeros(0, 2, zeros)

# Solving for the corresponding interpolating polynomial of degree n = 5
n = 5
points = zeros

# A matrix
A = A_matrix(n, points)

# b matrix
b = b_matrix(points, g)

# coefficients
coeffcheb_n5 = np.linalg.solve(A, b)

coeffcheb_n5_table = pd.DataFrame(coeffcheb_n5)
coeffcheb_n5_table.columns = ["coefficients"]
coeffcheb_n5_table.index = ['a5', 'a4', 'a3', 'a2', 'a1', 'a0']

# interpolating polynomial with points as the zeroes of the corresponding chebyshev polynomial
def polycheb_n5(x):
    return poly(5, coeffcheb_n5, x)



# FOR THE n = 10 INTERPOLATING POLYNOMIAL WITH POINTS EQUAL TO ZEROS OF m = n + 1 = 11 th degree ChebPoly:
# get the 11 th degree chebyshev polynomial
m = 11

# get the zeros of the cheb poly in [-1, 1]
zeros = cheb_zeros(m)

# transform zeros to interval [a, b] = [0, 2]
zeros = tfmd_zeros(0, 2, zeros)

# Getting the interpolating polynomial of degree n = 10 with points as zeros of 11th degree cheb poly
n = 10
points = zeros

# A matrix
A = A_matrix(n, points)

# b matrix
b = b_matrix(points, g)

# coefficients
coeffcheb_n10 = np.linalg.solve(A, b)

coeffcheb_n10_table = pd.DataFrame(coeffcheb_n10)
coeffcheb_n10_table.columns = ["coefficients"]
coeffcheb_n10_table.index = ['b10', 'b9', 'b8', 'b7', 'b6', 'b5', 'b4', 'b3', 'b2', 'b1', 'b0']

# Interpolating function
def polycheb_n10(x):
    return poly(10, coeffcheb_n10, x)


# // Plot the interpolants on the same graph //
x = np.linspace(0, 2, 100)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.spines['bottom'].set_position(('data', 0))
ax.spines['left'].set_position(('data', 0))
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

plt.plot(x, polycheb_n5(x), "b", label="cheb, n = 5")
plt.plot(x, polycheb_n10(x), "r", label="cheb, n = 10")
plt.plot(x, g(x))

plt.legend(loc="lower center")


# 3. // Compute the cubic splines interpolants using equally spaced interpolation points
# g(x) in [0, 2]
# FOR THE n = 5 CUBIC SPLINE INTERPOLANT (n+1 points)
n = 5
nr_points = n + 1
cubpoints_n5 = np.linspace(0, 2, nr_points)

interval = [0, 2]
h = (interval[1] - interval[0])/n

# Function for calculating A matrix for cubic splines problem Ax = b
# (solving for coefficients -the x vector: a0, ..., an- in the cubic splines linear combination)

def Acub_matrix(n):
    # first row
    zeroes = np.tile(0, n)
    r0 = np.append([1], zeroes)

    # middle rows + combine w/ first
    A = np.array([])
    for i in range(n - 1):
        lzeroes = np.tile(0, i)
        rzeroes = np.tile(0, (n + 1) - (len(lzeroes) + 3))
        new_row = np.append(np.append(lzeroes, [1, 4, 1]), rzeroes)
        if i == 0:
            A = np.append(A, new_row)
        else:
            A = np.vstack((A, new_row))

    A = np.vstack((r0, A))

    # last row + combine w/ previous
    rn = np.append(zeroes, [1])
    A = np.vstack((A, rn))
    return A


# Function for calculating the b matrix in cubic splines problem Ax = b
def bcub_matrix(p, h, func):
    coeff = [1] + [6]*(len(p) - 2) + [1]
    return [[coeff[i] / (h**3) * func(p[i])] for i in range(len(p))]



# A matrix
A = Acub_matrix(5)

# b matrix
b = bcub_matrix(cubpoints_n5, h, g)

# solving for a0, ..., an coefficients
coeff = np.linalg.solve(A, b)

# solving for a-1 and an+1 coefficients
a_l = 2*coeff[0] - coeff[1]
a_r = 2*coeff[n] - coeff[n - 1]

# vector of all coefficients from a-1, a0, ..., an, an+1
cubcoeff_n5 = np.vstack((a_l, np.vstack((coeff, a_r))))

cubcoeff_n5_table = pd.DataFrame(cubcoeff_n5)
cubcoeff_n5_table.columns = ["coefficients"]
cubcoeff_n5_table.index = ['a6', 'a5', 'a4', 'a3', 'a2', 'a1', 'a0', 'a-1']


# Defining the cubic splines base function B_0(x) - for the first point x0
# INPUTS:
#   - a = first point
#   - h = sub interval length
#   - x = variable
def B_0(x, a, h):
    return np.piecewise(x, [x <= a - 2*h, ((a - 2*h < x) & (x <= a - h)), ((a - h < x) & (x <= a)),
                            ((a < x) & (x <= a + h)), ((a + h < x) & (x <= a + 2*h)), x > a + 2*h],
                        [0, lambda x: 1/6 * (2*h + (x - a))**3,
                         lambda x: (2 * h**3)/3 - 0.5 * (x - a)**2 * (2*h + (x - a)),
                         lambda x: (2 * h**3)/3 - 0.5 * (x - a)**2 * (2*h - (x - a)),
                         lambda x: 1/6 * (2*h - (x - a))**3, 0])


# Defining general cubic spline base function B_k(x), k = -1, (0, 1, ..., n), n + 1
def B_k(x, a, h, k):
    return B_0(x - k*h, a, h)


# Calculating the complete cubic spline approximation for n = 5
def S3_n5(x, n=5, a=cubpoints_n5[0], h=0.4, coeff=cubcoeff_n5):
    S3 = 0
    for k in range(n + 3):
        S3 += coeff[k] * B_k(x, a, h, k - 1)
    return S3

# REPEATING FOR THE CUBIC SPLINE MADE UP OF SAMPLING n + 1 = 11 POINTS FOR n = 10 SUB INTERVALS
n = 10
nr_points = n + 1
cubpoints_n10 = np.linspace(0, 2, nr_points)


interval = [0, 2]
h = (interval[1] - interval[0]) / n

# A matrix
A = Acub_matrix(10)

# b matrix
b = bcub_matrix(cubpoints_n10, h, g)

# coefficients a0, ..., an
coeff = np.linalg.solve(A, b)

# computing a-1 and an+1 coefficients as well
a_l = 2*coeff[0] - coeff[1]
a_r = 2*coeff[n] - coeff[n - 1]

# vector of all coefficients in order
cubcoeff_n10 = np.vstack((a_l, np.vstack((coeff, a_r))))

cubcoeff_n10_table = pd.DataFrame(cubcoeff_n10)
cubcoeff_n10_table.columns = ["coefficients"]
cubcoeff_n10_table.index = ['b11', 'b10', 'b9', 'b8', 'b7', 'b6', 'b5', 'b4', 'b3', 'b2', 'b1', 'b0', 'b-1']

# cubic spline function for n = 10
def S3_n10(x, n=10, a=cubpoints_n10[0], h=0.2, coeff=cubcoeff_n10):
    S3 = 0
    for k in range(n + 3):
        S3 += coeff[k] * B_k(x, a, h, k - 1)
    return S3


# // Plot the interpolants on the same graph //
x = np.linspace(0, 2, 100)
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))

plt.plot(x, S3_n5(x), "b", label="cubic spline, n = 5")
plt.plot(x, S3_n10(x), "r", label="cubic spline, n = 10")
plt.plot(x, g(x))

plt.legend(loc="lower center")

# QUESTIONS TO ANSWER IN ESSAY:
# Briefly consider the differences between the two interpolants in each case
# What conclusions can you draw from these differences?

# differences: relative to original, shape, intersects w original, accuracy















# -------------------------------- Question 2: --------------------------------
# f(x) = (1 - x^2)sin(2*pi*x) on the interval [-1, 2]

# 1. // Compute the interpolating polynomials of degree n = 6 and n = 10 using equally spaced
#       interpolation points.
# GENERAL FUNCTION
def f(x):
    return (1 - x**2) * np.sin(2 * np.pi * x)

# COMPUTING INTERPOLATING POLYNOMIAL OF DEGREE n = 6
n = 6
nr_points = n + 1
polyn6points = np.linspace(-1, 2, nr_points)

# A matrix
A = A_matrix(n, polyn6points)

# b matrix
b = b_matrix(polyn6points, f)

# coefficients
coeff_n6 = np.linalg.solve(A, b)

# condition number
k = np.linalg.norm(np.linalg.inv(A)) * np.linalg.norm(A)
print(k)

# interpolating polynomial
def poly2_n6(x):
    return poly(6, coeff_n6, x)


# REPEATING FOR INTERPOLATING POLYNOMIAL OF DEGREE n = 10
n = 10
nr_points = n + 1
polyn10points = np.linspace(-1, 2, nr_points)

# A matrix
A = A_matrix(n, polyn10points)

# b matrix
b = b_matrix(polyn10points, f)

# coefficients
coeff2_n10 = np.linalg.solve(A, b)

coeff2_n10_table = pd.DataFrame(coeff2_n10)
coeff2_n10_table.index = ['b10', 'b9', 'b8', 'b7', 'b6', 'b5', 'b4', 'b3', 'b2', 'b1', 'b0']

# condition number
k = np.linalg.norm(np.linalg.inv(A)) * np.linalg.norm(A)
print(k)

def poly2_n10(x):
    return poly(10, coeff2_n10, x)



# // Plot the interpolants and f(x) on the same graph //
x = np.linspace(-1, 2, 100)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

plt.plot(x, poly2_n6(x), 'b', label="n = 6")
plt.plot(x, poly2_n10(x), 'r', label="n = 10")
plt.plot(x, f(x), 'black', label="original f(x)")

plt.legend(loc="lower right")


# 2. // Compute interpolating polynomials of degree n = 6 and n = 10 based on the zeros of
#       the appropriate Chebyshev polynomial //
# Procedure: Get zeros of m = n+1 th degree chebyshev polynomials, use as points for interpolating polynomial, calculate
#            its coefficients

# FOR n = 6 INTERPOLATING POLYNOMIAL BASED OFF ZEROS OF m = 7 th DEGREE CHEBYSHEV POLYNOMIAL:
m = 7

# zeros in interval [-1, 1]
zeros = cheb_zeros(m)
chebn6nodes = zeros

# transforming zeros to interval [a, b] = [-1, 2]
zeros = tfmd_zeros(-1, 2, zeros)

# FINDING THE n = 6th DEGREE INTERPOLATING POLYNOMIAL
n = 6
polychebn6points = zeros

# A matrix
A = A_matrix(n, polychebn6points)

# b matrix
b = b_matrix(polychebn6points, f)

# coefficients
coeffcheb_n6 = np.linalg.solve(A, b)

# condition number
k = np.linalg.norm(np.linalg.inv(A)) * np.linalg.norm(A)
print(k)


# Polynomial function
def polycheb2_n6(x):
    return poly(6, coeffcheb_n6, x)


# FOR THE n = 10 CASE, NEED m = 11th degree Cheb poly:
m = 11

# getting zeros in interval [-1, 1]
zeros = cheb_zeros(m)
chebn10nodes = zeros

# transforming zeros to [-1, 2]
zeros = tfmd_zeros(-1, 2, zeros)

# Calculating interpolating polynomial
n = 10
polychebn10points = zeros

# A matrix
A = A_matrix(n, polychebn10points)

# b matrix
b = b_matrix(polychebn10points, f)

# coefficients
coeffcheb2_n10 = np.linalg.solve(A, b)

# condition number
k = np.linalg.norm(np.linalg.inv(A)) * np.linalg.norm(A)
print(k)

# Polynomial function
def polycheb2_n10(x):
    return poly(10, coeffcheb2_n10, x)


# // Plot the interpolants and f(x) on the same graph //
x = np.linspace(-1, 2, 100)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))

plt.plot(x, polycheb2_n6(x), "b", label="cheb, n = 6")
plt.plot(x, polycheb2_n10(x), "r", label="cheb, n = 10")
plt.plot(x, f(x), "black", label="original f(x)")

plt.legend(loc="upper center")

# 3. // Compute the cubic splines interpolants using equally spaced interpolation points.
# FOR THE n = 6 CASE LOOK AT THE CUBIC SPLINE CREATED USING n + 1 = 7 EQUALLY SPACED POINTS
n = 6
nr_points = n + 1
cubpoints2_n6 = np.linspace(-1, 2, nr_points)

interval = [-1, 2]
h = (interval[1] - interval[0]) / n

# A matrix
A = Acub_matrix(6)

# b matrix
b = bcub_matrix(cubpoints2_n6, h, f)

# solving for coefficients of the points; a0, ..., an
coeff = np.linalg.solve(A, b)

# condition number
k = np.linalg.norm(np.linalg.inv(A)) * np.linalg.norm(A)
print(k)

# solving for leading/trailing coefficients; a-1, an+1
a_l = 2*coeff[0] - coeff[1]
a_r = 2*coeff[n] - coeff[n - 1]

# vector of all coefficients in order
cubcoeff2_n6 = np.vstack((a_l, np.vstack((coeff, a_r))))


# cubic spline for n = 6 (looking at n + 1 = 7 equally spaced points)
def S3_2n6(x, n=6, a=cubpoints2_n6[0], h=0.5, coeff=cubcoeff2_n6):
    S3 = 0
    for k in range(n + 3):
        S3 += coeff[k] * B_k(x, a, h, k - 1)
    return S3


# FOR THE n = 10 CASE: Cubic spline based on n + 1 = 11 equally spaced points
n = 10
nr_points = n + 1
cubpoints2_n10 = np.linspace(-1, 2, nr_points)

interval = [-1, 2]
h = (interval[1] - interval[0]) / n

# A matrix
A = Acub_matrix(10)

# b matrix
b = bcub_matrix(cubpoints2_n10, h, f)

# solving for coefficients of points; a0, ..., an
coeff = np.linalg.solve(A, b)

# condition number
k = np.linalg.norm(np.linalg.inv(A)) * np.linalg.norm(A)
print(k)

# solving for leading and trailing extra base spline coefficients; a-1, an+1
a_l = 2*coeff[0] - coeff[1]
a_r = 2*coeff[n] - coeff[n - 1]

# vector of all coefficients in cubic spline equation, in order
cubcoeff2_n10 = np.vstack((a_l, np.vstack((coeff, a_r))))

# Cubic spline approximation based on 11 equally spaced points on f
def S3_2n10(x, n=10, a=cubpoints2_n10[0], h=0.3, coeff=cubcoeff2_n10):
    S3 = 0
    for k in range(n + 3):
        S3 += coeff[k] * B_k(x, a, h, k - 1)
    return S3


# // Plot the interpolants and f(x) on the same graph //
x = np.linspace(-1, 2, 100)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))

plt.plot(x, S3_2n6(x), "b", label="cubic spline, n = 6")
plt.plot(x, S3_2n10(x), "r", label="cubic spline, n = 10")
plt.plot(x, f(x), "black", label="original f(x)")

plt.legend(loc="upper center")
plt.show()



# QUESTIONS TO ANSWER IN ESSAYS:
# 1. What is the actual maximum error of each interpolation? (calculate)
# 2. Compare the maximum error with the theoretical errors for interpolating polynomials and Chebyshev polynomials (calculate vs formulas)
# 3. Why are the interpolations using the n = 6 particularly bad in each case?
# 4. Explain briefly why an interpolating polynomial of degree less than seven will never provide a good approximation to this f(x)


# 1. Calculate maximum error as maximum magnitude of |function(x) - approx(x)|, i.e max(abs(func(x) - approx(x)))
# Max errors for n = 6 polynomial, cheb polynomial, and cubic spline approximations
maxerr_polyn6 = np.max(np.abs(f(x) - poly2_n6(x)))
maxerr_polychebn6 = np.max(np.abs(f(x) - polycheb2_n6(x)))
maxerr_cubn6 = np.max(np.abs(f(x) - S3_2n6(x)))

# Max errors for n = 10 polynomial, cheb polynomial, and cubic spline approximations
maxerr_polyn10 = np.max(np.abs(f(x) - poly2_n10(x)))
maxerr_polychebn10 = np.max(np.abs(f(x) - polycheb2_n10(x)))
maxerr_cubn10 = np.max(np.abs(f(x) - S3_2n10(x)))

print(maxerr_polyn6, maxerr_polychebn6, maxerr_cubn6, "\n", maxerr_polyn10, maxerr_polychebn10, maxerr_cubn10)

# 2. Theoretical errors for interpolating polynomials and Chebyshev polynomials (equations in notes)

def d7f(x):
    return 64*np.pi**5*(14*np.pi*x*np.sin(2*np.pi*x) + 2*np.pi**2*(x**2 - 1)*np.cos(2*np.pi*x) - 21*np.cos(2*np.pi*x))

def abs_d7f(x):
    return np.abs(64*np.pi**5*(14*np.pi*x*np.sin(2*np.pi*x) + 2*np.pi**2*(x**2 - 1)*np.cos(2*np.pi*x) - 21*np.cos(2*np.pi*x)))

def neg_abs_d7f(x):
    return -np.abs(64*np.pi**5*(14*np.pi*x*np.sin(2*np.pi*x) + 2*np.pi**2*(x**2 - 1)*np.cos(2*np.pi*x) - 21*np.cos(2*np.pi*x)))


# need n + 1 = 10 + 1 = 11th derivative of x
def d11f(x):
    return 1024*np.pi**9*(22*np.pi*x*np.sin(2*np.pi*x) + 2*np.pi**2*(x**2 - 1)*np.cos(2*np.pi*x) - 55*np.cos(2*np.pi*x))

def abs_d11f(x):
    return np.abs(1024*np.pi**9*(22*np.pi*x*np.sin(2*np.pi*x) + 2*np.pi**2*(x**2 - 1)*np.cos(2*np.pi*x) - 55*np.cos(2*np.pi*x)))

def neg_abs_d11f(x):
    return -np.abs(1024*np.pi**9*(22*np.pi*x*np.sin(2*np.pi*x) + 2*np.pi**2*(x**2 - 1)*np.cos(2*np.pi*x) - 55*np.cos(2*np.pi*x)))

# product in formula
def maxproduct(points, x = np.linspace(-1, 2, 100)):
    product = 1
    for i in range(len(points)):
        product = product * np.abs(x - points[i])
    return np.max(product)


# polynomial theoretical errors
n = 6
x_max = fminbound(neg_abs_d7f, 1.5, 2) # by inspecting the plot of abs_d7f below to search in restricted bound
theo_maxerr_polyn6 = maxproduct(polyn6points)/np.math.factorial(n + 1) * abs_d7f(x_max)
print(theo_maxerr_polyn6)

theo_maxerr_chebpolyn6 = (3/2)**(n+1) * 1/(2**n * np.math.factorial(n + 1)) * abs_d7f(x_max)
print(theo_maxerr_chebpolyn6)

n = 10
x_max = fminbound(neg_abs_d11f, 1.5, 2) # by inspecting the plot of abs_d7f below to search in restricted bound
theo_maxerr_polyn10 = maxproduct(polyn10points)/np.math.factorial(n + 1) * abs_d11f(x_max)
print(theo_maxerr_polyn10)

theo_maxerr_chebpolyn10 = (3/2)**(n+1) * 1/(2**n * np.math.factorial(n + 1)) * abs_d11f(x_max)
print(theo_maxerr_chebpolyn10)

# 3. The n = 6 interpolations are particularly bad because the point samples all have an f(x_k) value of 0, and so it is
#    approximated as a straight line y = 0.
# 4. Any cubic spline approximation based on more points will always be at least as accurate

x = np.linspace(-1, 2, 100)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))

# verifying that x_max points actually coincide with maximal absolute value
plt.plot(x, d7f(x), "b", label="abs d7f(x)")
plt.show()
plt.plot(x, d11f(x), "c", label="neg abs d11f(x)")
plt.show()

# Checking that the chebyshev interpolation for question 2, n = 10 is correct: error must be in the abs_d11f(x_max) part
actual = np.polynomial.chebyshev.Chebyshev.interpolate(lambda x: (1 - x**2)*np.sin(2*np.pi*x), 10, domain=[-1, 2])

plt.plot(x, actual(x))
plt.plot(x, f(x))
plt.show()