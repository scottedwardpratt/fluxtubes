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

a10=np.arange(0,11)
a20=np.arange(0,21)
a30=np.arange(0,31)
a40=np.arange(0,41)

mydata=np.loadtxt('../trajectories/A20.dat',skiprows=0,unpack=False)
for i in arange(0,20):
  plt.plot(a20,mydata[i],color='r',linestyle='-',linewidth=2,marker='None')

#amax=5
#z10=1.5*(amax-(a10-amax)*(a10-amax)/amax)
#plt.plot(a10,z10,color='k',linestyle='-',linewidth='6',marker='None')
amax=10
z20=1.5*(amax-(a20-amax)*(a20-amax)/amax)
plt.plot(a20,z20,color='k',linestyle='-',linewidth=5)
#amax=15
#z30=1.5*(amax-(a30-amax)*(a30-amax)/amax)
#plt.plot(a30,z30,color='k',linestyle='-',linewidth='6',marker='None')
#amax=20
#z40=1.5*(amax-(a40-amax)*(a40-amax)/amax)
#plt.plot(a40,z40,color='k',linestyle='-',linewidth='6',marker='None')

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,81,5), minor=False)
ax.set_xticklabels(np.arange(0,81,5), minor=False, family='serif')
ax.set_xticks(np.arange(0,81,1), minor=True)
#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.xaxis.set_major_formatter(sformatter)
plt.xlim(0.0,20)

ax.set_yticks(np.arange(0,200,20), minor=False)
ax.set_yticklabels(np.arange(0,200,20), minor=False, family='serif')
ax.set_yticks(np.arange(0,200,5), minor=True)
plt.ylim(0.0,50.0)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

text(15,13,'$A_{\\rm max}=20$')
text(30,24,'$A_{\\rm max}=40$')
#text(45,35,'$A_{\\rm max}=60$')
#text(60,46,'$A_{\\rm max}=80$')

plt.xlabel('$A$', fontsize=18, weight='normal')
plt.ylabel('$\langle C_1\\rangle$',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('traj.pdf',format='pdf')
os.system('open -a Preview traj.pdf')
#plt.show()
quit()
