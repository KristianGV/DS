import matplotlib.pyplot as plt
import numpy as np
from sympy.solvers import solve
from sympy import Symbol


def f1(x):
  return -np.sqrt(3)*x-1

def f2(x):
  return np.sqrt(3)*x-1

def f3(x):
  return 1/2
  
def S(x,y):
  lambda1=(1-2*x)/3
  lambda2=(x+np.sqrt(3)*y+1)/3
  lambda3=((x-np.sqrt(3)*y+1)/3)
  return -lambda1*np.log(lambda1)-lambda2*np.log(lambda2)-lambda3*np.log(lambda3)

theta=np.linspace(0,2*np.pi,num=200)
x_b=np.sin(theta)
y_b=np.cos(theta)

x = Symbol('x')
x1, =  solve(f1(x)-f2(x))
x2, =  solve(f1(x)-f3(x))
x3, =  solve(f2(x)-f3(x))
y1 = f1(x1)
y2 = f1(x2)
y3 = f2(x3)




plt.figure(1)

plt.fill([x1,x2,x3,x1],[y1,y2,y3,y1],'red',alpha=0.5,label='Positive eigenvalues',zorder=0)
plt.plot(x_b,y_b,color='black',linestyle='--',label='|m|=1')


plt.plot(x1,f1(x1))
plt.plot(x2,f1(x2))
plt.plot(x3,f2(x3))



m1=np.linspace(-1,1)
m8=np.linspace(-1,1)

X,Y=np.meshgrid(m1,m8)
Z=S(Y,X)
plt.contourf(X, Y, Z,alpha=1,zorder=5)
plt.colorbar()


plt.xlabel("$m_1$")
plt.ylabel("$m_8$")
plt.legend(loc=1,facecolor ='w',framealpha=1)
# plt.savefig('pos_eig_val.pdf',format='pdf',dpi=1000)
plt.show()
