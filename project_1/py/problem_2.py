import numpy as np
import matplotlib.pyplot as plt
import math

## variables using Earth and the Sun for now
m1 = 10.0**24
m2 = 10.0**30
rm = (m1*m2)/(m1+m2)
r = 100000000.0
theta = 0.0
rdot = 1.0
thetadot = 0.0
p_r = rm*rdot
p_theta = rm*(r**2)*thetadot

G = 6.67408*10**(-11)
k = G*m1*m2

## step size
a = 0.0     # Start of the interval
b = 10.0    # End of the interval
N = 1000 # Number of steps
h = (b-a)/N # Size of a single step


## create arrays
t_points = np.arange(a,b,h)
p_r_points = []
p_theta_points = []
r_points = []
theta_points = []

## iterate through the algorithm
for t in t_points:
    p_theta = p_theta
    theta = theta + h*((p_theta)/(rm*(r**2)))
    p_r = p_r - h*(k/(r**2))
    r = r + (p_r)/rm

    p_r_points.append(p_r)
    r_points.append(r)
    p_theta_points.append(p_theta)
    theta_points.append(theta)



## plot
plt.plot(t_points, r_points, "b-")
plt.xlabel("t")
plt.ylabel("r")
plt.title("Attempt at numerical position plot")
plt.show()
