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
c3c1=mydata[15]
c3c1error=mydata[16]
c4c3=mydata[17]
c4c3error=mydata[18]

def parabola(y, a, b):
    return a*y**2 + b
def gaussian(y, A, sigma):
    return A*np.exp(-y**2/(2*sigma**2))/(np.sqrt(2*np.pi*abs(sigma)))
def quartic(y, A, B, C):
    return A*y**4+B*y**2+C
def hexic(y, A, B, C, D):
    return A*y**6+B*y**4+C*y**2+D

fig=plt.figure(figsize=(8,12))
ax1 = fig.add_subplot(211)
plt.ylabel('$C_3/C_1$', fontsize=22, weight='normal')
ax2 = fig.add_subplot(212)
plt.ylabel('$C_4/C_3$', fontsize=22, weight='normal')


param, param_cov = curve_fit(parabola, y, c3c1)
parab = param[0]*(y)**2+param[1]
param, param_cov = curve_fit(quartic, y, c3c1)
quart = param[0]*y**4+param[1]*y**2+param[2]
#param, param_cov = curve_fit(hexic, y, c3c1)
#hex = param[0]*y**6+param[1]*y**4+param[2]*y**2+param[3]
ax1.errorbar(y,c3c1,c3c1error,linestyle='',linewidth=2,markersize=8,color='b', marker='d', markerfacecolor=None, markeredgecolor=None)
ax1.plot(y,parab,linestyle='-',color='b')
ax1.plot(y,quart,linestyle=':',color='b')
#ax1.plot(y,hex,linestyle='--',color='b')

param, param_cov = curve_fit(parabola, y, c4c3)
parab = (param[0]*(y)**2+param[1])
param, param_cov = curve_fit(quartic, y, c4c3)
quart = param[0]*y**4+param[1]*y**2+param[2]
ax2.errorbar(y,c4c3,c4c3error,linestyle='',linewidth=2,markersize=8,color='k', marker='^', markerfacecolor=None, markeredgecolor=None)
ax2.plot(y,parab,linestyle='-',color='k')
ax2.plot(y,quart,linestyle=':',color='k')
plt.xlabel('y',fontsize=18 , weight='normal')

plt.savefig('figs/altratios.pdf',format='pdf')
os.system('xdg-open figs/altratios.pdf')
quit()
