import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

file="edist.dat";
mydata = np.loadtxt(file,skiprows=1,unpack=True)

edens=mydata[0]
n=mydata[1]

def gaussian(y, A, B, sigma):
    return A*np.exp(-(y-B)**2/(2*sigma**2))/(np.sqrt(2*np.pi*abs(sigma)))
def fit1(y, A, sigma):
    return A*y**2*np.exp(-y**2/(2*sigma**2))/(np.sqrt(2*np.pi*abs(sigma)))
def fit2(y, A, r):
    return A*y**4*np.exp(-y/r)

fig=plt.figure(figsize=(8,12))
param, param_cov = curve_fit(gaussian, edens, n, p0=[1.,25.,10.])
gauss=param[0]*np.exp(-(edens-param[1])**2/(2*param[2]**2))/(np.sqrt(2*np.pi*abs(param[2])))
print(param)
param, param_cov = curve_fit(fit1, edens, n, p0=[1.,30.])
eq1=param[0]*edens**2*np.exp(-edens**2/(2*param[1]**2))/(np.sqrt(2*np.pi*abs(param[1])))
print(param)
param, param_cov = curve_fit(fit2, edens, n, p0=[10.,20.])
eq2=param[0]*edens**4*np.exp(-edens/param[1])
print(param)
plt.plot(edens,n,linestyle='',linewidth=2,markersize=8,color='k', marker='d', markerfacecolor=None, markeredgecolor=None)
plt.plot(edens,eq1,linestyle='-',color='r',label='$x^2$*gaussian')
plt.plot(edens,eq2,linestyle='-',color='g',label='$x^4 e^{-x}$')
plt.plot(edens,gauss,linestyle='-',color='b',label='gaussian')
plt.xlabel('dE/dy',fontsize=18 , weight='normal')
plt.ylabel('% of paths',fontsize=18 , weight='normal')
plt.legend(fontsize=18)

plt.savefig('figs/dist.pdf',format='pdf')
os.system('xdg-open figs/dist.pdf')
quit()
