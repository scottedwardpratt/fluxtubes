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
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,6))
fig = plt.figure(1)
x0=0.145
y0=0.1
plotwidth=1.0-x0-0.03
plotheight=(0.98-y0)/3.0


B0=1
Bf=1
Ng=6

#############################################
ipanel=0
itraj=14
filename='trajectories/B0'+str(B0)+'Bf'+str(Bf)+'Ng'+str(Ng)+'/traj'+str(itraj)+'.txt'
print('filename=',filename)
mydata = np.loadtxt(filename,skiprows=0,unpack=True)
#mydata = np.loadtxt('trajectories/B01Bf1Ng6/traj0.txt',skiprows=0,unpack=True)
ytilde=mydata[0]
#print(ytilde)
ytilde=-1+2.0*(ytilde)/(Ng+1)
print(ytilde)
Q2=mydata[3]
pminusq=mydata[1]-mydata[2]
Q3=mydata[4]

ax = fig.add_axes([x0,y0+ipanel*plotheight,plotwidth,plotheight])

plt.plot(ytilde,Q2,lw=4,drawstyle='steps',color='g',label="$Q^{(2)}_{\\rm left}$")
plt.plot(ytilde,pminusq,lw=4,drawstyle='steps',ls='--',color='r',label="$(p-q)_{\\rm left}$")

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(-2,2,0.5), minor=False)
ax.set_xticklabels(np.arange(-2,2,0.5), minor=False, family='serif')
ax.set_xticks(np.arange(-2,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
plt.xlim(-1.0,1.0)

ax.set_yticks(np.arange(-50,50,5), minor=False)
ax.set_yticklabels(np.arange(-50,50,5), minor=False, family='serif')
ax.set_yticks(np.arange(-50,50,1), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
plt.ylim(-8,15)

plt.xlabel('$\\tilde{y}$', fontsize=18, weight='normal')
plt.ylabel(None)
 
#ax.legend(framealpha=100, loc="lower left",fontsize='18')
#############################################
#############################################
ipanel=1
itraj=15
filename='trajectories/B0'+str(B0)+'Bf'+str(Bf)+'Ng'+str(Ng)+'/traj'+str(itraj)+'.txt'
print('filename=',filename)
mydata = np.loadtxt(filename,skiprows=0,unpack=True)
#mydata = np.loadtxt('trajectories/B01Bf1Ng6/traj0.txt',skiprows=0,unpack=True)
ytilde=mydata[0]
#print(ytilde)
ytilde=-1+2.0*(ytilde)/(Ng+1)
print(ytilde)
Q2=mydata[3]
pminusq=mydata[1]-mydata[2]
Q3=mydata[4]

ax = fig.add_axes([x0,y0+ipanel*plotheight,plotwidth,plotheight])

plt.plot(ytilde,Q2,lw=4,drawstyle='steps',color='g',label="$Q^{(2)}_{\\rm left}$")
plt.plot(ytilde,pminusq,lw=4,drawstyle='steps',ls='--',color='r',label="$(p-q)_{\\rm left}$")

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(-2,2,0.5), minor=False)
#ax.set_xticklabels(np.arange(-2,2,0.5), minor=False, family='serif')
ax.set_xticklabels([])
ax.set_xticks(np.arange(-2,2,0.1), minor=True)
plt.xlim(-1.0,1.0)

ax.set_yticks(np.arange(-50,50,5), minor=False)
ax.set_yticklabels(np.arange(-50,50,5), minor=False, family='serif')
ax.set_yticks(np.arange(-50,50,1), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
plt.ylim(-8,15)

plt.xlabel(None)
plt.ylabel('$(p-q)_{\\rm left}$, $Q^{(2)}_{\\rm left}$',fontsize=18)
 
#ax.legend(framealpha=100, loc="lower left",fontsize='18')
#############################################
#############################################
ipanel=2
itraj=17
filename='trajectories/B0'+str(B0)+'Bf'+str(Bf)+'Ng'+str(Ng)+'/traj'+str(itraj)+'.txt'
print('filename=',filename)
mydata = np.loadtxt(filename,skiprows=0,unpack=True)
#mydata = np.loadtxt('trajectories/B01Bf1Ng6/traj0.txt',skiprows=0,unpack=True)
ytilde=mydata[0]
#print(ytilde)
ytilde=-1+2.0*(ytilde)/(Ng+1)
print(ytilde)
Q2=mydata[3]
pminusq=mydata[1]-mydata[2]
Q3=mydata[4]

ax = fig.add_axes([x0,y0+ipanel*plotheight,plotwidth,plotheight])

plt.plot(ytilde,Q2,lw=4,drawstyle='steps',color='g',label="$Q^{(2)}_{\\rm left}$")
plt.plot(ytilde,pminusq,lw=4,drawstyle='steps',ls='--',color='r',label="$(p-q)_{\\rm left}$")

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xticks(np.arange(-2,2,0.5), minor=False)
#ax.set_xticklabels(np.arange(-2,2,0.5), minor=False, family='serif')
ax.set_xticklabels([])
ax.set_xticks(np.arange(-2,2,0.1), minor=True)
plt.xlim(-1.0,1.0)

ax.set_yticks(np.arange(-50,50,5), minor=False)
ax.set_yticklabels(np.arange(-50,50,5), minor=False, family='serif')
ax.set_yticks(np.arange(-50,50,1), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
plt.ylim(-8,15)

plt.xlabel(None)
plt.ylabel(None)
 
ax.legend(loc="upper right",fontsize='16',framealpha=0.5,ncol=2)
#############################################

plt.savefig('sample_trajectories.pdf',format='pdf')
os.system('open -a Preview sample_trajectories.pdf')
