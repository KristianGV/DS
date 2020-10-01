import matplotlib.pyplot as plt
import numpy as np

#dataset= np.genfromtxt('examples/aux/generic_wimp_oh2-planck-sigmav.dat')
M,oh2 =np.loadtxt('generic_FIMP.dat',unpack=True,skiprows=1, usecols=range(0,2)) 

plt.figure(1)
plt.semilogx(10/M,oh2, color='red',linestyle='dotted' )
plt.xlabel('$m_\chi/m_{Z_p}$')
plt.ylabel('$\Omega h^2$')
plt.show()