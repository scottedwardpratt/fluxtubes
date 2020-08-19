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
Ksigma=mydata[13]
Kerror=mydata[14]

def parabola(y, a):
    return a*y**2 - 49*a
def gaussian(y, A, sigma):
    return A*np.exp(-y**2/(2*sigma**2))/(np.sqrt(2*np.pi*abs(sigma)))
def quartic(y, A, B):
    return A*y**4+B*y**2 - (49*49*A+49*B)
def hexic(y, A, B, C):
    return A*y**6+B*y**4+C*y**2 - (A*7**6+B*7**4+C*7**2)
def octic(y, A, B, C, D):
    return A*y**8+B*y**6+C*y**4+D*y**2 - (A*7**8+B*7**6+C*7**4+D*7**2)

fig=plt.figure(figsize=(8,16))
ax1 = fig.add_subplot(411)
plt.ylabel('$C_1$', fontsize=22, weight='normal')
ax2 = fig.add_subplot(412)
plt.ylabel('$C_2$', fontsize=22, weight='normal')
ax3 = fig.add_subplot(413)
plt.ylabel('$C_3$', fontsize=22, weight='normal')
ax4 = fig.add_subplot(414)
plt.ylabel('$C_4$', fontsize=22, weight='normal')


param, param_cov = curve_fit(parabola, y, kappa1_avg)
parab = param[0]*(y**2-49)
ax1.errorbar(y,kappa1_avg,k1error,linestyle='',linewidth=2,markersize=8,color='r', marker='o', markerfacecolor=None, markeredgecolor=None)
ax1.plot(y,parab,linestyle='-',color='r')

param, param_cov = curve_fit(quartic, y, kappa2_avg)
quart = param[0]*y**4+param[1]*y**2 - (49*49*param[0]+49*param[1])
param, param_cov = curve_fit(parabola, y, kappa2_avg)
parab = param[0]*(y**2-49)
ax2.errorbar(y,kappa2_avg,k2error,linestyle='',linewidth=2,markersize=8,color='g', marker='d', markerfacecolor=None, markeredgecolor=None)
ax2.plot(y,quart,linestyle=':',color='g')

param, param_cov = curve_fit(hexic, y, kappa3_avg)
hex = param[0]*y**6+param[1]*y**4+param[2]*y**2 - (param[0]*7**6+param[1]*7**4+param[2]*7**2)
param, param_cov = curve_fit(quartic, y, kappa3_avg)
quart = param[0]*y**4+param[1]*y**2 - (49*49*param[0]+49*param[1])
ax3.errorbar(y,kappa3_avg,k3error,linestyle='',linewidth=2,markersize=10,color='b', marker='^', markerfacecolor=None, markeredgecolor=None)
ax3.plot(y,quart,linestyle=':',color='b')
ax3.plot(y,hex,linestyle='--',color='b')

param, param_cov = curve_fit(quartic, y, kappa4_avg)
quart = param[0]*y**4+param[1]*y**2 - (49*49*param[0]+49*param[1])
param, param_cov = curve_fit(hexic, y, kappa4_avg)
hex = param[0]*y**6+param[1]*y**4+param[2]*y**2 - (param[0]*7**6+param[1]*7**4+param[2]*7**2)
param, param_cov = curve_fit(octic, y, kappa4_avg)
oct = param[0]*y**8+param[1]*y**6+param[2]*y**4+param[3]*y**2 - (param[0]*7**8+param[1]*7**6+param[2]*7**4+param[3]*7**2)
ax4.errorbar(y,kappa4_avg,k4error,linestyle='',linewidth=2,markersize=10,color='k', marker='s', markerfacecolor=None, markeredgecolor=None)
ax4.plot(y,quart,linestyle=':',color='k')
ax4.plot(y,hex,linestyle='--',color='k')
ax4.plot(y,oct,linestyle='-',color='k')
plt.xlabel('y',fontsize=18 , weight='normal')

plt.savefig('figs/moments.pdf',format='pdf')
os.system('xdg-open figs/moments.pdf')
quit()
