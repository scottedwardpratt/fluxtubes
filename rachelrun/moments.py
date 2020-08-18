import matplotlib.pyplot as plt
import numpy as np
import os

file="moments.dat";
mydata = np.loadtxt(file,skiprows=1,unpack=True)

Amax=mydata[0]
Ssigma=mydata[1]
Serror=mydata[2]
Ksigma=mydata[3]
Kerror=mydata[4]
omega=mydata[5]
werror=mydata[6]


fig=plt.figure(figsize=(8,12))
ax1 = fig.add_subplot(311)
plt.ylabel('$\\omega=C_2/C_1$', fontsize=22, weight='normal')
ax2 = fig.add_subplot(312)
plt.ylabel('$S\sigma=C_3/C_2$', fontsize=22, weight='normal')
ax3 = fig.add_subplot(313)
plt.ylabel('$K\sigma^2=C_4/C_2$', fontsize=22, weight='normal')
ax1.errorbar(Amax,omega,werror,linestyle='-',linewidth=2,markersize=8,color='r', marker='o', markerfacecolor=None, markeredgecolor=None)
ax2.errorbar(Amax,Ssigma,Serror,linestyle='-',linewidth=2,markersize=8,color='g', marker='^', markerfacecolor=None, markeredgecolor=None)
ax3.errorbar(Amax,Ksigma,Kerror,linestyle='-',linewidth=2,markersize=10,color='b', marker='s', markerfacecolor=None, markeredgecolor=None)
plt.xlabel('A',fontsize=18 , weight='normal')

plt.savefig('moments.pdf',format='pdf')
os.system('xdg-open moments.pdf')
quit()
