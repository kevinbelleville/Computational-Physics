import numpy as np
import matplotlib.pyplot as plt

nt = 30
temp = np.linspace(1,4, nt)
crit_temp = 2.26

def onsager(t):
    if t < crit_temp:
        return (1-(np.sinh(2/t))**(-4))**(1/8)
    else:
        return 0

o = []
for _t in temp:
    o.append(onsager(_t))

plt.plot(temp, o)
plt.show()
