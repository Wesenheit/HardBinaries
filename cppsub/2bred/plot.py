import matplotlib.pyplot as plt
import numpy as np

A=np.loadtxt("file.txt")
plt.plot(A[:,0],A[:,1])
plt.savefig("plot.png")
