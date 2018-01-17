import numpy as np
import matplotlib.pyplot as plt
from values import *

name = "SE2"

def f(x):
    return m*(omega**2)*x

def g(p):
    return p/m



t_points = np.arange(a,b,h)
x_points = []
v_points = []
e_points = []

def se2(x,p):
    for t in t_points:
        e = 0.5*(1/m)*((p**2)+(omega**2)*(x**2))
        e_points.append(e)
        x_points.append(x)
        v_points.append(p)
        x = x + h*g(p)
        p = p - h*f(x)

"""
plt.plot(t_points, x_points, "r-")
plt.plot(x_points, p_points, "g-")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.show()
"""
