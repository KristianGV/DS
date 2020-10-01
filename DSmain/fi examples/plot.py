import matplotlib.pyplot as plt
import numpy as np

#dataset= np.genfromtxt('examples/aux/generic_wimp_oh2-planck-sigmav.dat')
tmp,dydt =np.loadtxt('dydt.dat',unpack=True,skiprows=1, usecols=range(0,2)) 

plt.figure(1)
plt.semilogx(tmp,dydt, color='red')
plt.xlabel('$T$ [GeV]')
plt.ylabel('$dY/d\log_{10}T$')
plt.show()