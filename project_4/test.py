import matplotlib.pyplot as plt 
import numpy as np 

a = [1,2,3,4,5,4,3,2,1]
t = [0,1,2,3,4,5,6,7,8]

corr = np.correlate(a,a, "full")


plt.plot(t, a)
plt.show()

plt.plot([i for i in range(len(corr))], corr)
plt.show()

# testing
# test = np.correlate(sunspots_normalized, sunspots_normalized, "full")
# test = test[test.size//2:]
# test_x = [i for i in range(len(test))]
# plt.plot(test_x, test)
# plt.show()