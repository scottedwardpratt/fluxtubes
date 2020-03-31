import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

mydata10 = np.loadtxt('Cbar_A10.dat',skiprows=0,unpack=True)
a10=mydata10[0]
cbar10=mydata10[1]
mydata20 = np.loadtxt('Cbar_A20.dat',skiprows=0,unpack=True)
a20=mydata20[0]
cbar20=mydata20[1]
mydata30 = np.loadtxt('Cbar_A30.dat',skiprows=0,unpack=True)
a30=mydata30[0]
cbar30=mydata30[1]
mydata40 = np.loadtxt('Cbar_A40.dat',skiprows=0,unpack=True)
a40=mydata40[0]
cbar40=mydata40[1]



plt.plot(a10,cbar10,color='r',linestyle='None',markersize=4, marker='o', markerfacecolor='r', markeredgecolor='r')
plt.plot(a20,cbar20,color='r',linestyle='None',markersize=4, marker='o', markerfacecolor='r', markeredgecolor='r')
plt.plot(a30,cbar30,color='r',linestyle='None',markersize=4, marker='o', markerfacecolor='r', markeredgecolor='r')
plt.plot(a40,cbar40,color='r',linestyle='None',markersize=4, marker='o', markerfacecolor='r', markeredgecolor='r')

amax=10
z10=1.5*(amax-(a10-amax)*(a10-amax)/amax)
plt.plot(a10,z10,color='k',linestyle='-',linewidth='6',marker='None')
amax=20
z20=1.5*(amax-(a20-amax)*(a20-amax)/amax)
plt.plot(a20,z20,color='k',linestyle='-',linewidth='6',marker='None')
amax=30
z30=1.5*(amax-(a30-amax)*(a30-amax)/amax)
plt.plot(a30,z30,color='k',linestyle='-',linewidth='6',marker='None')
amax=40
z40=1.5*(amax-(a40-amax)*(a40-amax)/amax)
plt.plot(a40,z40,color='k',linestyle='-',linewidth='6',marker='None')

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,81,20), minor=False)
ax.set_xticklabels(np.arange(0,81,10), minor=False, family='serif')
ax.set_xticks(np.arange(0,81,1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,80)

ax.set_yticks(np.arange(0,200,20), minor=False)
ax.set_yticklabels(np.arange(0,200,20), minor=False, family='serif')
ax.set_yticks(np.arange(0,200,5), minor=True)
plt.ylim(0.0,64.0)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

text(15,13,'$A_{\\rm max}=20$')
text(30,24,'$A_{\\rm max}=40$')
text(45,35,'$A_{\\rm max}=60$')
text(60,46,'$A_{\\rm max}=80$')

plt.xlabel('$A$', fontsize=18, weight='normal')
plt.ylabel('$\langle C_1\\rangle$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('cbar.pdf',format='pdf')
os.system('open -a Preview cbar.pdf')
#plt.show()
quit()
