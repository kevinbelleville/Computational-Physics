## Import all of my methods
import euler
import rk2
import rk4
import se1
import se2
import sho1

import numpy as np
import matplotlib.pyplot as plt

## Execute the calculations of each one
euler.euler(euler.x, euler.v)
rk2.rk2(rk2.x, rk2.v)
rk4.rk4(rk4.x, rk4.v)
se1.se1(se1.x, se1.p)
se2.se2(se2.x, se2.p)

methods = [euler, rk2, rk4, se1, se2]
colors = ["r-", "y-", "g-", "b-", "c-"]

## Graphing area
for i, method in enumerate(methods):
    plt.plot(method.t_points, method.e_points, colors[i], label = method.name)
    plt.legend()
    plt.xlabel("t")
    plt.ylabel("total energy")
    plt.title("Energy Conservation")

plt.plot(sho1.t_points, sho1.e_points, "r--", label = sho1.name)
plt.show()
