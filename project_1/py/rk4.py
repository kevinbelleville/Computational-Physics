import numpy as np
import matplotlib.pyplot as plt
from values import *

name = "RK4"

def f(x):
    return -omega**2 * x


t_points = np.arange(a,b,h)
x_points = []
v_points = []
e_points = []

def rk4(x,v):
    for t in t_points:
        e = 0.5*m*((v**2)+(omega**2)*(x**2))
        e_points.append(e)
        x_points.append(x)
        v_points.append(v)
        v1 = h*f(x)
        x1 = h*v
        v2 = h*f(x+0.5*v1)
        v3 = h*f(x+0.5*v2)
        v4 = h*f(x+v3)
        x += x1
        v += (v1+2*v2+2*v3+v4)/6

"""
plt.plot(t_points, x_points, "r-")
plt.plot(x_points, v_points, "g-")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.show()
"""
