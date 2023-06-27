# Testing Jacobi, Gauss-Seidel & Successive relaxation methods.
import numpy as np
import math as m

# -------------------------- Code for question 1 of the assignment ---------------------------------
print("----------------------- Question 1 -----------------------")


# // Create initial problem elements Ax = b -> A, b, D, L, U, initial estimate x, and iteration matrices //
# Function that creates a vector with elements args written as "x1, x2, x3"
def vector(*args):
    return np.array([[x] for x in args])


# Defining all preliminary parts
A = np.array([[9, 1, 2, 1, 0, 4], [0, 8, 1, 2, 1, 0], [1, 0, 9, 1, 2, 1], [0, 1, 0, 7, 1, 2], [0, 0, 1, 0, 9, 1],
              [4, 0, 0, 1, 0, 6]])

b = vector(27, 22, 17, 23, 13, 24)

D = np.diag(np.diag(A))

U = -np.triu(A, k=1)

L = -np.tril(A, k=-1)

# Defining iteration matrices for each method: Jacobi, Gauss-Seidel, Successive Relaxation
B_J = np.matmul(np.linalg.inv(D), L + U)
B_GS = np.matmul(np.linalg.inv(D - L), U)


def B_SR(w):
    return np.matmul(np.linalg.inv(D - w * L), (1 - w) * D + w * U)


# Defining function giving spectral radius of input matrix x
def spectral_radius(B):
    return np.max(np.abs(np.linalg.eigvals(B)))


# // i) determining optimal w value for the matrix //
# Given that B_SR is a 6x6 matrix with non-zero diagonal, w_optimal will be in range (0,2) -theorem-
# calculate spectral radii for w values in this range
w = np.linspace(0, 2, 1000)
v = []
for i in range(len(w) - 1):
    v.append(spectral_radius(B_SR(w[i])))

# find the optimal w value
w_optimal = w[np.argmin(v)]
print(f"The optimal w value is: \n {w_optimal}")

# // ii) Finding spectral radii for the Jacobi, Gauss-Seidel and Successive Relaxation methods
J_specrad = spectral_radius(B_J)
print(f"\nThe spectral radius for the Jacobi method is: \n {J_specrad}")
GS_specrad = spectral_radius(B_GS)
print(f"The spectral radius for the Gauss-Seidel method is: \n {GS_specrad}")
SR_specrad = spectral_radius(B_SR(w_optimal))
print(f"The spectral radius for the Successive Relaxation method is: \n {SR_specrad}")

# // iii) Number of iterations required to obtain a solution of accuracy at least 10^-8 for J, GS, SR methods //
# exact solution of x to Ax = b
x = np.linalg.solve(A, b)

# actual solution is x = (1, 2, 1, 2, 1, 3) so an appropriate guess would be:
x_initial_guess = vector(1, 1, 1, 1, 1, 1)


# Jacobi iterative method
def J_iterations(x, x_guess):
    v = [x_guess]
    for i in range(100):
        v.append(np.matmul(B_J, v[i]) + np.matmul(np.linalg.inv(D), b))
        error = np.linalg.norm(x - v[i+1], np.inf)
        if error <= 10 ** -8:
            print(f"\nThe number of iterations for the Jacobi method is: \n {len(v) - 1}")
            return (len(v) - 1)
            break


# Gauss-Seidel iterative method
def GS_iterations(x, x_guess):
    v = [x_guess]
    for i in range(100):
        v.append(np.matmul(B_GS, v[i]) + np.matmul(np.linalg.inv(D - L), b))
        error = np.linalg.norm(x - v[i], np.inf)
        if error <= 10 ** -8:
            print(f"The number of iterations for the Gauss-Seidel method is: \n {len(v) - 1}")
            return (len(v) - 1)
            break


# Successive Relaxation iterative method
def SR_iterations(x, x_guess, w):
    v = [x_guess]
    for i in range(100):
        v.append(np.matmul(B_SR(w), v[i]) + np.matmul(np.linalg.inv(D - w * L), w * b))
        error = np.linalg.norm(x - v[i], np.inf)
        if error <= 10 ** -8:
            print(f"The number of iterations for the Successive Relaxation method is: \n {len(v) - 1}")
            return (len(v) - 1)
            break


J_iterations(x, x_initial_guess)
GS_iterations(x, x_initial_guess)
SR_iterations(x, x_initial_guess, w_optimal)

# --------------------------------- Code for question 2 of assignment -----------------------------------
print("\n----------------------- Question 2 -----------------------")
# // Create initial problem elements Ax = b -> A, b, D, L, U, initial estimate x, and iteration matrices //
# Defining all preliminary parts
A = np.array([[7, 2, 0, 0, 0, 0], [5, 4, 3, 0, 0, 0], [0, 3, 9, 1, 0, 0], [0, 0, 3, 7, 1, 0], [0, 0, 0, 3, 9, 1],
              [0, 0, 0, 0, 5, 6]])

b = vector(20, 34, 47, 27, 17, 17)

D = np.diag(np.diag(A))

U = -np.triu(A, k=1)

L = -np.tril(A, k=-1)

# Redefining the iteration matrices
B_J = np.matmul(np.linalg.inv(D), L + U)


def B_SR(w):
    return np.matmul(np.linalg.inv(D - w * L), (1 - w) * D + w * U)


# // i) determining optimal w value for the matrix //
# Given that B_SR is a 6x6 matrix with non-zero diagonal, w_optimal will be in range (0,2) -theorem-
# calculate spectral radii for w values in this range
w = np.linspace(0, 2, 1000)
v = []
for i in range(len(w) - 1):
    v.append(spectral_radius(B_SR(w[i])))

# calculate optimal w
w_optimal = w[np.argmin(v)]
print(f"The optimal w value is: \n {w_optimal}")

# // ii) finding spectral radii for Jacobi and Successive Relaxation methods, how much faster does SR converge? //
J_specrad = spectral_radius(B_J)
print(f"\nThe spectral radius for the Jacobi method is: \n {J_specrad}")
SR_specrad = spectral_radius(B_SR(w_optimal))
print(f"The spectral radius for the Successive Relaxation method is: \n {SR_specrad}")

# calculating how much faster SR method is
r = m.log(SR_specrad, J_specrad)
print(f"\nThe Successive Relaxation method is expected to converge {r} times faster than the Jacobi method")

# // iii) Number of iterations required to obtain a solution of accuracy at least 10^-8 for J, GS, SR methods //
# exact solution of x to Ax = b
x = np.linalg.solve(A, b)

# actual solution is x = (1, 3, 3, 2, 1, 2) so will use the same approximate x_initial_guess = (1, 1, 1, 1, 1, 1)

# counting number of iterations
J_it = J_iterations(x, x_initial_guess)
SR_it = SR_iterations(x, x_initial_guess, w_optimal)

# see how much faster the SR method converges
print(f"\nThe Successive Relaxation method actually converges {J_it / SR_it} times faster")
