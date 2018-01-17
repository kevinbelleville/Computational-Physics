import math
import matplotlib.pyplot as plt
import numpy as np
from values import *

name = "Analytical Solution"

def f(t, omega, x0):
    return x0*math.cos(omega*t)

def g(t, omega, x0):
    return -x0*omega*math.sin(omega*t)

x0 = 1.0
x1 = 2.0
x2 = 5.0
x3 = 0.5

x_points = []
v_points = []
t_points = np.arange(a,b,h)

e_points = []

for t in t_points:
    e = 0.5*m*((v**2)+(omega**2)*(x**2))
    e_points.append(e)
    x_points.append(f(t, omega, x0))
    v_points.append(g(t, omega, x0))

x1points = []
v1points = []

for t in t_points:
    x1points.append(f(t, omega, x1))
    v1points.append(g(t, omega, x1))

x2points = []
v2points = []

for t in t_points:
    x2points.append(f(t, omega, x2))
    v2points.append(g(t, omega, x2))

x3points = []
v3points = []

for t in t_points:
    x3points.append(f(t, omega, x3))
    v3points.append(g(t, omega, x3))

stuff2 = [x_points, x1points, x2points, x3points]
stuff = [v_points, v1points, v2points, v3points]
colors = ["r-", "b-", "c-", "g-"]
names = ["x = 1.0", "x = 2.0", "x = 5.0", "x = 0.5"]



def sho1(stuff = stuff, stuff2 = stuff2):
    for i in range(4):
        plt.plot(stuff[i], stuff2[i], colors[i], label = names[i])
    plt.legend()
    plt.title("Fixed Frequency, Variable Initial x_0")
    plt.xlabel('v')
    plt.ylabel('x')
    plt.show()
