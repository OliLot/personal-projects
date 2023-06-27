# Part 1
from numpy import *
from pylab import *

# Defining the function
def f(x):
  return 0.01 - 0.4*x + (x*x)/(1+x*x)

# Plotting the function to find positive roots

# Range of x values covering roots
x = arange(0,2.5,0.01)

# The plot
plot(x,f(x))
xlabel('x')
ylabel('f(x)')
title('f(x) on limited range')
# Displaying the figure.
plt.grid(True)
show()

# Thus, know there are three roots: One close to 0.1, another close to 0.4,
# and a final one close to 2

# Newton-Raphson scheme
def dfdx(x):
  return -0.4 + (2*x*(1+x*x)-x*x*2*x)/(1+x*x)**2

#Solving scheme for all three roots
def NRroot(x,f,dfdx):
    for i in range(10):
        x = x - f(x)/dfdx(x)
        z.append(round(x,4))
        print(round(x,4))
        if z[i+1] == z[i]:
            print(" ")
            print("Iteration n =",i,"yields the correct result as confirmed by iteration n =",i+1)
            print("The difference; xn - xn-1 =",z[i]-z[i-1],"where xn-1 =",z[i-1])
            break
    return x

xguess = [0.1, 0.4, 2]
for i in range(3):
    print("Applying NR method and solving for positive root number", i+1,":")
    z = [0]
    r = NRroot(xguess[i],f,dfdx)
    print("and the positive root xn =",round(r,4))
    print(" ")

#Moving on to next question
print('')
print('--------------------------------------------------------------------')
print('')

# RK method
def f(x,t):
  return 0.01 - 0.4*x + (x*x)/(1+x*x)

# Fourth order Runge-Kutta scheme to integrate function for single time step:
def RK4step(x,t,h,f):
  k1 = h*f(x,t)
  k2 = h*f(x+k1/2.0,t+h/2.0)
  k3 = h*f(x+k2/2.0,t+h/2.0)
  k4 = h*f(x+k3,t+h)
  return x+(k1+2*k2+2*k3+k4)/6

# Choice of timestep, h, and total number of steps to take
h = 0.5
nsteps = 80

# a) Initial condition x(0) = 0.45
x = 0.45
t = 0.0

# Creating arrays in which to store solutions
xstore = []
tstore = []
xstore.append(x)
tstore.append(t)

# Now solve the ODE using the Runge-Kutta method
for n in range(nsteps+1):
  x = RK4step(x,t,h,f)
  t = t+h
  xstore.append(x)
  tstore.append(t)
  
# Plotting the results  
plot(tstore,xstore,'k-')
xlabel('t')
ylabel('x')
ylim(0,0.5)
title('x(t) for x(0)=0.45')
show()

# Note: same step size and number of steps are used (a & b). They result in a
# plot accurately depicting x vs t (which is in the range of 0 < t < 40),
# because for t>40 both plots have a horizontal asymptote. Picking a smaller
# step size and a corresponding greater number of steps (to fit the range of t
# values) yields a plot that looks identical

# b) Initial condition = x(0)=0.5
x = 0.5
t = 0.0

# Creating arrays in which to store solutions
xstore = []
tstore = []
xstore.append(x)
tstore.append(t)

# Solving the ODE using the Runge-Kutta method
for n in range(nsteps+1):
  x = RK4step(x,t,h,f)
  t = t+h
  xstore.append(x)
  tstore.append(t)
  
# Plotting the results  
plot(tstore,xstore,'k-')
xlabel('t')
ylabel('x')
ylim(0,2.2)
title('x(t) for x(0)=0.50')
show()

#--------------------------------------------------------------------------------------------------

# Part 2

from pylab import *
import numpy as np

# Question 1

# Have coupled system:
# 1-x-0.1*exp(0.1*y)=0 and
# 1-y-0.2*exp(0.3*x)=0, so plotting both and seeing where they are equal:

def f(x):
    return np.log((1-x)/0.1)/0.1

def g(x):
    return 1-0.2*exp(0.3*x)

# range up to 1 as x < 1 is the domain for f(x)
x = arange(0,1,0.01)

# The plot
plot(x,f(x))
plot(x,g(x))
xlabel('x')
ylabel('y')
title('f(x) and g(x) on limited range')
# Displaying the figure.
plt.grid(True)
show()

# Thus, a rough guess for the coordinates of the fourth steady state is
# (0.83, 1)

# Define this as a vector valued function:
def f(x,y):
    return array([[1-x-0.1*exp(0.1*y)],[1-y-0.2*exp(0.3*x)]])

# 2D Newton-Raphson scheme to find steady states
def J(x,y):
    return array([[-1,-0.01*exp(0.1*y)],[-0.06*exp(0.3*x),-1]])

def MultiNRroot(x,y,f,J):
    for i in range(10):
        [[x],[y]]=[[x],[y]]-np.matmul(np.linalg.inv(J(x,y)),f(x,y))
        z = [[x],[y]]
        print("x =",str(z)[2:-24], "y=",str(z)[24:-2])
    return [[x],[y]]

print("Performing 10 iterations of the 2D NR-scheme with initial guess (0.83,1):")
r = MultiNRroot(0.83,1,f,J)
print(" ")
print("Thus, the fourth steady state, to an accuracy of 16 decimal places is at")
print("x =",str(r[0])[1:-1],"and y =",str(r[1])[1:-1])


print('---------------------------------------')

# Question 2

# Now take the same RHS to integrate the system
def f(x,y):
    return array([x*(1-x-0.1*exp(0.1*y)),y*(1-y-0.2*exp(0.3*x))])

# Fourth order RK scheme for one time step
def RK4step(x,y,t,h,f):
  k1 = h*f(x,y)
  k2 = h*f(x+k1[0]/2.0,y+k1[1]/2.0)
  k3 = h*f(x+k2[0]/2.0,y+k2[1]/2.0)
  k4 = h*f(x+k3[0],y+k3[1])
  return [x+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,y+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6]

#Define the time step and number of steps to use
h=0.1
nsteps=1000

# initial conditions
x_initial = [0, 0.0001, 0.1, 0.4, 0.7, 0.892, 1.1, 1.5, 2, 5, 0.05, 0.05, 0.05, 0.25, 0.892, 2, 2, 2]
y_initial = [0, 0.8, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0.25, 0.05, 0.05, 0.05, 0.15, 0.5, 0.739]
t = 0

# Solving the coupled ODEs for this initial condition
for i in range(18):
    xstore = []
    ystore = []
    xstore.append(x_initial[i])
    ystore.append(y_initial[i])
    for n in range(nsteps+1):
        [x,y] = RK4step(xstore[n],ystore[n],t,h,f)
        t = t+h
        xstore.append(x)
        ystore.append(y)
    plot(xstore, ystore, 'k-')
    xlabel('x')
    ylabel('y')
    title('Phase plane plot')
    xlim(0,2)
    ylim(0,2)
        
#Showing steady states (0,0), (0, 0.8), (0.9, 0), (0.892, 0.739)
plot(0,0,color='red', marker='o', markersize=7)
plot(0,0.8,color='red', marker='o', markersize=7)
plot(0.9,0,color='red', marker='o', markersize=7)
plot(0.892,0.739,color='red', marker='o', markersize=7)

print("Note: The steady states are shown as red dots on the produced phase plane plot. The (0,0.8) steady state is an unstable saddlepoint (with trajectories approaching along the y axis, and diverging for x>0). The steady state (0,0) is an unstable steady point. Steady point (0.9,0) is also a saddle point (has trajectories approaching along the x axis, but diverging for y>0). The steady point (0.892, 0.739) is stable.") 