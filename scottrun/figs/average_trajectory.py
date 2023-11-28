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
        'size'   : 18}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(15,4.5))
fig = plt.figure(1)
x0=0.07
plotwidth=(1.0-x0-0.01)/3.0

plt.rcParams['ytick.right'] = True
plt.rcParams['ytick.labelright'] = False
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.labelleft'] = True

#############################################
#############################################

B0=1
Bf=1
Ng=12
filename='results/B0'+str(B0)+'Bf'+str(Bf)+'Ng'+str(Ng)+'.txt'
print('filename=',filename)
mydata = np.loadtxt(filename,skiprows=0,unpack=True)
#mydata = np.loadtxt('trajectories/B01Bf1Ng6/traj0.txt',skiprows=0,unpack=True)
ytilde=mydata[0]
print(ytilde)
ytilde=-1+2.0*(ytilde-3*B0+1)/(Ng+1)
print(ytilde)
Q2=mydata[1]
pminusq=mydata[3]
Q3=mydata[2]
#############################################

ax = fig.add_axes([x0,0.16,plotwidth,0.8])

plt.plot(ytilde,Q2,drawstyle='steps',linestyle='-',lw=4,color='g',label="$\langle Q^{(2)}\\rangle_{\\rm left}$")
plt.plot(ytilde,10*pminusq,drawstyle='steps',linestyle='--',lw=4,color='r',label="$\langle p-q\\rangle_{\\rm left}\\times 10$")


besta=0.25
bestgmag=1.66667

#print('ytilde=',ytilde)
#bestguess=100000000.0
#for a in arange(0.2,0.3,0.001):
#  for gmag in arange(1.7,2.5,0.001):
#    guess=-gmag*a*sinh(ytilde/a)
#    offby=0.0
#    for i in range(3,15):
#      offby=offby+(guess[i]-10*pminusq[i])*(guess[i]-10*pminusq[i])
    
#    if offby < bestguess:
#      bestguess=offby
#      print ('besta=',besta,' gmag=',bestgmag,' offby=',offby)

#      besta=a
#      bestgmag=gmag

guess=-bestgmag*besta*sinh(ytilde/besta)   
plt.plot(ytilde,guess,lw=3,color='b')
#offbyvec=(guess-10*pminusq)
#offby=np.linalg.norm(offbyvec)
#print ('offby=',offby)
#print(guess)
#print(10*pminusq)

#plt.plot(B,4*B,linestyle='--',lw=1,color='b',label="$\langle Q^{(2)}\\rangle$")

#plt.plot(B,Q3,linestyle='dotted',lw=4,color='r',label="$\langle Q^{(3)}\\rangle$")
#plt.plot(B,pminusq,linestyle='-',lw=4,color='g',label="$\langle p-q\\rangle$")

ax.tick_params(axis='both', which='major', labelsize=18, direction='inout')
ax.set_xticks(np.arange(-2,2,0.4), minor=False)
ax.set_xticklabels(np.arange(-2,2,0.4), minor=False, family='serif')
ax.set_xticks(np.arange(-2,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
plt.xlim(-1,1)

ax.set_yticks(np.arange(-50,50,5), minor=False)
ax.set_yticklabels(np.arange(-50,50,5), minor=False, family='serif')
ax.set_yticks(np.arange(-50,50,1), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0d'))
plt.ylim(-10,15)

plt.xlabel('$\\tilde{y}$', fontsize=22, weight='normal')
plt.ylabel('10\\times$\langle p-q\\rangle_{\\rm left}$, $\langle Q^{(2)}\\rangle_{\\rm left}$',fontsize=24,labelpad=-3)

ax.legend(framealpha=100, loc="lower left",fontsize='20')

bigtitle='$B_{\\rm target}=$'+str(B0)+', $B_{\\rm proj}=$'+str(Bf)+', $N_g$='+str(Ng);
littletitle='howdy'
text(0,13,bigtitle,fontsize='22',ha='center')




#############################################
#second panel
#############################################

B0=1
Bf=1
Ng=6
filename='results/B0'+str(B0)+'Bf'+str(Bf)+'Ng'+str(Ng)+'.txt'
print('filename=',filename)
mydata = np.loadtxt(filename,skiprows=0,unpack=True)
#mydata = np.loadtxt('trajectories/B01Bf1Ng6/traj0.txt',skiprows=0,unpack=True)
ytilde=mydata[0]
print(ytilde)
ytilde=-1+2.0*(ytilde-3*B0+1)/(Ng+1)
print(ytilde)
Q2=mydata[1]
pminusq=mydata[3]
Q3=mydata[2]
#############################################

ax = fig.add_axes([x0+plotwidth,0.16,plotwidth,0.8])

plt.plot(ytilde,Q2,drawstyle='steps',linestyle='-',lw=4,color='g',label="$\langle Q^{(2)}\\rangle_{\\rm left}$")
plt.plot(ytilde,10*pminusq,drawstyle='steps',linestyle='--',lw=4,color='r',label="$\langle p-q\\rangle_{\\rm left}\\times 10$")

#plt.plot(B,4*B,linestyle='--',lw=1,color='b',label="$\langle Q^{(2)}\\rangle$")

#plt.plot(B,Q3,linestyle='dotted',lw=4,color='r',label="$\langle Q^{(3)}\\rangle$")
#plt.plot(B,pminusq,linestyle='-',lw=4,color='g',label="$\langle p-q\\rangle$")

ax.tick_params(axis='both', which='major', labelsize=18, direction='inout')
ax.set_xticks(np.arange(-2,2,0.4), minor=False)
ax.set_xticklabels(np.arange(-2,2,0.4), minor=False, family='serif')
ax.set_xticks(np.arange(-2,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
plt.xlim(-1,1)

ax.set_yticks(np.arange(-50,50,5), minor=False)
ax.set_yticklabels([])
ax.set_yticks(np.arange(-50,50,1), minor=True)
plt.ylim(-10,15)

plt.xlabel('$\\tilde{y}$', fontsize=22, weight='normal')
#plt.ylabel('$\langle p-q\\rangle_{\\rm left}$, $\langle Q^{(2)}\\rangle_{\\rm left}$',fontsize=18,labelpad=-3)
plt.ylabel(None)


bigtitle='$B_{\\rm target}=$'+str(B0)+', $B_{\\rm proj}=$'+str(Bf)+', $N_g$='+str(Ng);
littletitle='howdy'
text(0,13,bigtitle,fontsize='22',ha='center')


#############################################
# 3rd panel
#############################################

B0=4
Bf=1
Ng=6
filename='results/B0'+str(B0)+'Bf'+str(Bf)+'Ng'+str(Ng)+'.txt'
print('filename=',filename)
mydata = np.loadtxt(filename,skiprows=0,unpack=True)
#mydata = np.loadtxt('trajectories/B01Bf1Ng6/traj0.txt',skiprows=0,unpack=True)
ytilde=mydata[0]
print(ytilde)
ytilde=-1+2.0*(ytilde-3*B0+1)/(Ng+1)
print(ytilde)
Q2=mydata[1]
pminusq=mydata[3]
Q3=mydata[2]
#############################################

ax = fig.add_axes([x0+2*plotwidth,0.16,plotwidth,0.8])

plt.plot(ytilde,Q2,drawstyle='steps',linestyle='-',lw=4,color='g',label="$\langle Q^{(2)}\\rangle_{\\rm left}$")
plt.plot(ytilde,10*pminusq,drawstyle='steps',linestyle='--',lw=4,color='r',label="$\langle p-q\\rangle_{\\rm left}\\times 10$")

#plt.plot(B,4*B,linestyle='--',lw=1,color='b',label="$\langle Q^{(2)}\\rangle$")

#plt.plot(B,Q3,linestyle='dotted',lw=4,color='r',label="$\langle Q^{(3)}\\rangle$")
#plt.plot(B,pminusq,linestyle='-',lw=4,color='g',label="$\langle p-q\\rangle$")

ax.tick_params(axis='both', which='major', labelsize=18, direction='inout')
ax.set_xticks(np.arange(-2,2,0.4), minor=False)
ax.set_xticklabels(np.arange(-2,2,0.4), minor=False, family='serif')
ax.set_xticks(np.arange(-2,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1.1f'))
plt.xlim(-1,1)

ax.set_yticks(np.arange(-50,50,5), minor=False)
ax.set_yticklabels([])
ax.set_yticks(np.arange(-50,50,1), minor=True)
plt.ylim(-10,15)

plt.xlabel('$\\tilde{y}$', fontsize=22, weight='normal')
#plt.ylabel('$\langle p-q\\rangle_{\\rm left}$, $\langle Q^{(2)}\\rangle_{\\rm left}$',fontsize=18,labelpad=-3)
plt.ylabel(None)

bigtitle='$B_{\\rm target}=$'+str(B0)+', $B_{\\rm proj}=$'+str(Bf)+', $N_g$='+str(Ng);
littletitle='howdy'
text(0,13,bigtitle,fontsize='22',ha='center')


#############################################
plt.savefig('average_trajectory.pdf',format='pdf')
os.system('open -a Preview average_trajectory.pdf')
