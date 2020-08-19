import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

file="moments.dat";
mydata = np.loadtxt(file,skiprows=1,unpack=True)

y=mydata[0]
kappa1_avg=mydata[1]
k1error=mydata[2]
kappa2_avg=mydata[3]
k2error=mydata[4]
kappa3_avg=mydata[5]
k3error=mydata[6]
kappa4_avg=mydata[7]
k4error=mydata[8]
omega=mydata[9]
werror=mydata[10]
Ssigma=mydata[11]
Serror=mydata[12]
Ksigma2=mydata[13]
Kerror=mydata[14]

def parabola(y, a, b):
    return a*y**2 + b
def gaussian(y, A, sigma):
    return A*np.exp(-y**2/(2*sigma**2))/(np.sqrt(2*np.pi*abs(sigma)))
def quartic(y, A, B, C):
    return A*y**4+B*y**2+C

fig=plt.figure(figsize=(8,12))
ax1 = fig.add_subplot(311)
plt.ylabel('$C_2/C_1$', fontsize=22, weight='normal')
ax2 = fig.add_subplot(312)
plt.ylabel('$C_3/C_2$', fontsize=22, weight='normal')
ax3 = fig.add_subplot(313)
plt.ylabel('$C_4/C_2$', fontsize=22, weight='normal')


param, param_cov = curve_fit(parabola, y, omega)
parab = param[0]*(y)**2+param[1]
ax1.errorbar(y,omega,werror,linestyle='',linewidth=2,markersize=8,color='g', marker='d', markerfacecolor=None, markeredgecolor=None)
ax1.plot(y,parab,linestyle='-',color='g')

param, param_cov = curve_fit(parabola, y, Ssigma)
parab = (param[0]*(y)**2+param[1])
ax2.errorbar(y,Ssigma,Serror,linestyle='',linewidth=2,markersize=8,color='b', marker='^', markerfacecolor=None, markeredgecolor=None)
ax2.plot(y,parab,linestyle='-',color='b')

param, param_cov = curve_fit(parabola, y, Ksigma2)
parab = (param[0]*(y)**2+param[1])
param, param_cov = curve_fit(quartic, y, Ksigma2)
quart = param[0]*y**4+param[1]*y**2+param[2]
ax3.errorbar(y,Ksigma2,Kerror,linestyle='',linewidth=2,markersize=10,color='k', marker='s', markerfacecolor=None, markeredgecolor=None)
ax3.plot(y,parab,linestyle='-',color='k')
ax3.plot(y,quart,linestyle=':',color='k')
plt.xlabel('y',fontsize=18 , weight='normal')

plt.savefig('figs/ratios.pdf',format='pdf')
os.system('xdg-open figs/ratios.pdf')
quit()
