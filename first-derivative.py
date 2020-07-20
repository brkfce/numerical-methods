# a numerical method to calculate the first derivative of a function
# with comparison to an analytical solution
# this works best as a jupyter notebook, but this is just for personal use



# import libraries
import numpy as np
from math import *
import matplotlib.pyplot as plt



# initialise a space-dependent sin function

# initial parameters
xmax = 10                       # physical domain (m)
nx = 200                        # number of samples
dx = xmax / (nx - 1)            # grid increment dx (m)
x = np.linspace(0, xmax, nx)    # space coordinates

# initialise sin funcntion
l = 20 * dx                     # wavelength
k = 2 * pi / l                  # wavenumber
f = np.sin(k*x)                 # function



# plot sin function
plt.plot(x, f)
plt.title("Sin function")
plt.xlabel("x, m")
plt.ylabel("Amplitude")
plt.xlim(0, xmax)
plt.grid()
plt.show()



# first derivative with two points

# initialisation of numerical and analytical derivatives
nder = np.zeros(nx)             # numerical derivative
ader = np.zeros(nx)             # analytical derivative

# numerical derivative of a given function
for i in range(1, nx - 1):
    nder[i] = (f[i+1] - f[i-1]) / (2 * dx)

# anaytical derivative of a given function
ader = k * np.cos(k * x)
# exclude boundaries
ader[0] = 0
ader[nx -1] = 0

# error in nder (rms)
rms = np.sqrt(np.mean(nder - ader)**2)



# plot the derivative error
plt.plot(x, nder, label = "Numerical Derivative, 2 points", marker = "+", color = "blue")
plt.plot(x, ader, label = "Analytical Derivative", lw = 2, ls = "-", color = "black")
plt.plot(x, nder - ader, label = "Difference", lw = 2, ls = ":")
plt.title("First Derivative, Err (rms) - %.6f " % (rms) )
plt.xlabel("x, m")
plt.ylabel("Amplitude")
plt.legend(loc = "lower left")
plt.grid()
plt.show()



# plotting the number of points per wavelength
plt.plot(x, nder, label = "Numerical Derivative, 2 points", marker = "+", color = "blue")
plt.title("First Derivative, Error = %.6f, $n_\lambda$ = %.2f" % (rms, l/dx))
plt.xlabel("x, m")
plt.ylabel("Amplitude")
plt.legend(loc = "lower left")
plt.xlim((xmax / 2 - l, xmax / 2 + l))
plt.grid()
plt.show()



# investigate error as a function of grid points per wavelength

# define a range of number of points per wavelength
# loop over points, calculate corresponding wavelength and calculate error

# initialise vectors
nmin = 3
nmax = 16
na = np.zeros(nmax - nmin + 1)  #vector with number of points per wavelength
err = np.zeros(nmax - nmin + 1) #vector with error

j = -1                              # array index

# loop through finite difference derivative calculation
for n in range(nmin, nmax):

    j = j + 1
    na[j] = n

    # initialise sin function
    l = na[j] * dx              # wavelength    
    k = 2 * pi / l              # wavenumber
    f = np.sin(k * x)           # sin function

    # numerical derivative of the sin function
    for i in range(1, nx - 1):
        nder[i] = (f[i+1] - f[i-1]) / (2 * dx)

    # analytical detivative of the sin function
    ader = k * np.cos(k * x)
    # exclude boundaries
    ader[0] = 0
    ader[nx - 1] = 0

    i0 = np.int(nx / 2)
    # error (rms)
    err[j] = (nder[i0] - ader[i0])**2 / ader[i0]**2 * 100

# plotting error as a function of number of points per wavelength
plt.plot(na, err, ls = "-", color = "blue")
plt.title("Error as a function of $n_\lambda$ ")
plt.xlabel("n$_\lambda$")
plt.ylabel("rms")
plt.grid()
plt.show()
