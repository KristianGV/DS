import matplotlib.pyplot as plt
import numpy as np

#dataset= np.genfromtxt('examples/aux/generic_wimp_oh2-planck-sigmav.dat')
x,y =np.loadtxt('2to2rhs.dat',unpack=True,skiprows=1, usecols=range(0,2)) 

plt.figure(1)
plt.plot(x,y, color='red')
plt.xlabel('x')
plt.ylabel('RHS')
plt.show()