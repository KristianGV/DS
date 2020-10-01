import matplotlib.pyplot as plt
import numpy as np

#dataset= np.genfromtxt('examples/aux/generic_wimp_oh2-planck-sigmav.dat')

# mdm_dec,oh2_dec =np.loadtxt('scalar_FIMP_decay.dat',unpack=True,skiprows=1, usecols=range(0,2)) 

mdm_2to2,oh2_2to2 =np.loadtxt('scalar_FIMP_2to2.dat',unpack=True,skiprows=1, usecols=range(0,2)) 

mdm_mo, oh2_mo=np.loadtxt('micrOMEGA_data.dat',unpack=True,skiprows=1, usecols=range(0,2), delimiter=';') 

# x,y=np.loadtxt('debug.dat',unpack=True,skiprows=1, usecols=range(0,2)) 

lgmass=np.linspace(np.log(0.1e-2),np.log(0.3e6),140)
mass=np.exp(lgmass)


# oh2_dec_interp=np.interp(mass,mdm_dec,oh2_dec)
# oh2_mo_interp=np.interp(mass,mdm_mo,oh2_mo)
# to2to=oh2_mo_interp-oh2_interp


plt.figure(1)

##### DECAY
# plt.loglog(mdm_dec,oh2_dec, color='black',linestyle='--',label='my data:  decay') 

#### 2->2
plt.loglog(mdm_2to2,oh2_2to2, color='black',linestyle='--',label='my data:  2->2') 

#### MICROMEGAS DATA
# plt.loglog(mdm_mo,oh2_mo, color='red',linestyle='--',label='micrOMEGAs data') 
# plt.loglog(mass,to2to, color='blue',label='2->2 data') 


plt.xlabel('$M_{DM}$')
plt.ylabel('$\Omega h^2$')
plt.title('Freeze-in abundance')
plt.legend()
# plt.savefig('fig/scalar_fimp.pdf',format='pdf')
plt.show()

#####  DEBUG  ######
# plt.figure(2)
# plt.plot(x,y,color='black',linestyle='--')
# plt.title('Debug')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

######  DECAY ABUNDANCE PER TEMPERATURE  ###########
# T,Y =np.loadtxt('scalar_FIMP_Y.dat',unpack=True,skiprows=1, usecols=range(0,2)) 
# plt.figure(2)
# plt.loglog(T,Y,color='black', linestyle='--')
# plt.xlim(1000000,1)
# plt.xlabel('T')
# plt.ylabel('Y')
# plt.show()