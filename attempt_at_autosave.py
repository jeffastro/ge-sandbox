from matplotlib import pyplot as plt
import numpy as np

a = np.linspace(0,10)
b = np.linspace(20,30)

plt.figure()
plt.scatter(a,b)
plt.savefig("../052416/attempt.png")