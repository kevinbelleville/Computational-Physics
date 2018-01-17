import numpy as np
import matplotlib.pyplot as plt
import sho1
from values import *

name = "Euler's Method"

def f(x):
    return -omega**2 * x


t_points = np.arange(a,b,h)
x_points = []
v_points = []
e_points = []

def euler(x, v):
    for t in t_points:
        e = 0.5*m*((v**2)+(omega**2)*(x**2))
        e_points.append(e)
        x_points.append(x)
        v_points.append(v)

        v1 = h*f(x)
        x1 = h*v
        x += x1
        v += v1

"""
plt.plot(t_points, x_points, "g-")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.show()
"""
