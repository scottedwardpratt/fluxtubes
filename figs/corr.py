import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
#import scipy.special as sp
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
ax = fig.add_axes([0.18,0.12,0.75,0.8])

mydata10 = np.loadtxt('corr_A10.dat',skiprows=0,unpack=True)
a10=mydata10[0]
corr10=mydata10[1]
mydata20 = np.loadtxt('corr_A20.dat',skiprows=0,unpack=True)
a20=mydata20[0]
corr20=mydata20[1]
mydata30 = np.loadtxt('corr_A30.dat',skiprows=0,unpack=True)
a30=mydata30[0]
corr30=mydata30[1]
mydata40 = np.loadtxt('corr_A40.dat',skiprows=0,unpack=True)
a40=mydata40[0]
corr40=mydata40[1]


amax=20
plt.plot(a10/amax,corr10,color='r',linestyle='None',markersize=4, marker='o', markerfacecolor='r', markeredgecolor='r')
amax=40
plt.plot(a20/amax,corr20,color='g',linestyle='None',markersize=4, marker='o', markerfacecolor='g', markeredgecolor='g')
amax=60
plt.plot(a30/amax,corr30,color='b',linestyle='None',markersize=4, marker='o', markerfacecolor='b', markeredgecolor='b')
amax=80
plt.plot(a40/amax,corr40,color='cyan',linestyle='None',markersize=4, marker='o', markerfacecolor='cyan', markeredgecolor='cyan')

amax=80
z40=1.0+0.26*exp(-4.5*a40/amax)
plt.plot(a40/amax,z40,color='k',linestyle='-',linewidth='6',marker='None')

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,1.2,0.2), minor=False)
ax.set_xticklabels(np.arange(0,1.2,0.2), minor=False, family='serif')
ax.set_xticks(np.arange(0,1.2,0.05), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,1.0)

ax.set_yticks(np.arange(0,2.0,0.1), minor=False)
ax.set_yticklabels(np.arange(0,2.0,0.1), minor=False, family='serif')
ax.set_yticks(np.arange(0,2.0,0.02), minor=True)
plt.ylim(0.95,1.3)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\delta A/A_{\\rm max}$', fontsize=18, weight='normal')
plt.ylabel('$\\frac{\langle C_1(A_{\\rm max}/2-\delta A/2)C_1(A_{\\rm max}/2+\delta A/2)\\rangle}{\langle C_1(A_{\\rm max}/2-\delta A/2)\\rangle\langle C_1(A_{\\rm max}/2+\delta A/2)\\rangle}$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('corr.pdf',format='pdf')
os.system('open -a Preview corr.pdf')
#plt.show()
quit()
