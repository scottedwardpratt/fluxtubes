import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 15}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(5.25,4.5))
fig = plt.figure(1)
x0=0.16
plotwidth=(1.0-x0-0.01)

plt.rcParams['ytick.right'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.labelleft'] = True

#############################################
#############################################

B0=1
Bf=1
Ng=12
#filename='results/B0'+str(B0)+'Bf'+str(Bf)+'Ng'+str(Ng)+'.txt'
filename='bcorr.txt'
mydata = np.loadtxt(filename,skiprows=1,unpack=True)
delNg=mydata[0]
Bcorr=mydata[1]

xguess=arange(1,14)
yguess=0.72*0.562**xguess
#############################################

ax = fig.add_axes([x0,0.16,plotwidth,0.81])

plt.plot(delNg,Bcorr,linestyle='None',marker='o',color='g')
plt.plot(xguess,yguess,linestyle='-',marker='None',color='g')


#plt.plot(B,4*B,linestyle='--',lw=1,color='b',label="$\langle Q^{(2)}\\rangle$")

#plt.plot(B,Q3,linestyle='dotted',lw=4,color='r',label="$\langle Q^{(3)}\\rangle$")
#plt.plot(B,pminusq,linestyle='-',lw=4,color='g',label="$\langle p-q\\rangle$")

ax.tick_params(axis='both', which='major', labelsize=15, direction='inout')
ax.set_xticks(np.arange(0.0,20.0,4.0), minor=False)
ax.set_xticklabels(np.arange(0.0,20.0,4.0), minor=False, family='serif')
ax.set_xticks(np.arange(0.0,20.0,1.0), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
plt.xlim(0,13)

ax.set_yticks(np.arange(-1.2,1.2,0.3), minor=False)
ax.set_yticklabels(np.arange(-1.2,1.2,0.3))
ax.set_yticks(np.arange(-1.2,1.2,0.1), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
plt.ylim(-0.05,0.7)

plt.xlabel('$\Delta N_g$', fontsize=22, weight='normal')
#plt.ylabel('$\langle p-q\\rangle_{\\rm left}$, $\langle Q^{(2)}\\rangle_{\\rm left}$',fontsize=18,labelpad=-3)
plt.ylabel('$B(\Delta N_g)$',fontsize=22)

#bigtitle='$B_{\\rm target}=$'+str(B0)+', $B_{\\rm proj}=$'+str(Bf)+', $N_g$='+str(Ng);
#text(0,13,bigtitle,fontsize='22',ha='center')


#############################################
plt.savefig('bcorr.pdf',format='pdf')
os.system('open -a Preview bcorr.pdf')
