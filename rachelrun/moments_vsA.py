import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

file="moments_vsA.dat";
mydata = np.loadtxt(file,skiprows=1,unpack=True)

A=mydata[0]
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
Ksigma=mydata[13]
Kerror=mydata[14]

fig=plt.figure(figsize=(8,16))
ax1 = fig.add_subplot(411)
plt.ylabel('$C_1$', fontsize=22, weight='normal')
ax2 = fig.add_subplot(412)
plt.ylabel('$C_2$', fontsize=22, weight='normal')
#ax3 = fig.add_subplot(313)
#plt.ylabel('$C_2/C_1$', fontsize=22, weight='normal')
ax3 = fig.add_subplot(413)
plt.ylabel('$C_3$', fontsize=22, weight='normal')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax4 = fig.add_subplot(414)
plt.ylabel('$C_4$', fontsize=22, weight='normal')

plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlabel('A',fontsize=18 , weight='normal')

def linear(A, m):
    return m*A
def parabola(A, m):
    return m*A**2
def cubic(A, m):
    return m*A**3
def quartic(A, m):
    return m*A**4

param, param_cov = curve_fit(linear, A, kappa1_avg)
line = param[0]*A
ax1.errorbar(A,kappa1_avg,k1error,linestyle='',linewidth=2,markersize=8,color='r', marker='o', markerfacecolor=None, markeredgecolor=None)
ax1.plot(A,line,linestyle='-',color='r')

param, param_cov = curve_fit(parabola, A, kappa2_avg)
parab = param[0]*A**2
ax2.errorbar(A,kappa2_avg,k2error,linestyle='',linewidth=2,markersize=8,color='g', marker='d', markerfacecolor=None, markeredgecolor=None)
ax2.plot(A,parab,linestyle='--',color='g')

"""
param, param_cov = curve_fit(linear, A, omega)
line = param[0]*A
ax3.errorbar(A,omega,werror,linestyle='',linewidth=2,markersize=10,color='b', marker='^', markerfacecolor=None, markeredgecolor=None)
ax3.plot(A,line,linestyle='-',color='b')
"""

param, param_cov = curve_fit(parabola, A, kappa3_avg)
parab = param[0]*A**2
param, param_cov = curve_fit(cubic, A, kappa3_avg)
cub = param[0]*A**3
ax3.errorbar(A,kappa3_avg,k3error,linestyle='',linewidth=2,markersize=10,color='b', marker='^', markerfacecolor=None, markeredgecolor=None)
ax3.plot(A,parab,linestyle='--',color='b')
ax3.plot(A,cub,linestyle=':',color='b')

param, param_cov = curve_fit(parabola, A, kappa4_avg)
parab = param[0]*A**2
param, param_cov = curve_fit(cubic, A, kappa4_avg)
cub = param[0]*A**3
param, param_cov = curve_fit(quartic, A, kappa4_avg)
quart = param[0]*A**4
ax4.errorbar(A,kappa4_avg,k4error,linestyle='',linewidth=2,markersize=10,color='k', marker='s', markerfacecolor=None, markeredgecolor=None)
ax4.plot(A,parab,linestyle='--',color='k')
ax4.plot(A,cub,linestyle=':',color='k')
ax4.plot(A,quart,linestyle='-',color='k')

plt.savefig('figs/moments_vsA.pdf',format='pdf')
os.system('xdg-open figs/moments_vsA.pdf')
quit()
