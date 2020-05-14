import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fits
import sys
from scipy.optimize import curve_fit

#find a way to do this using the parameters generated in the other script later
A = int(sys.argv[1])
traj_total = np.loadtxt("/home/thomas/fluxtubes/thomasrun/trajectories/A" + str(A) +".txt")
avg = np.average(traj_total, axis=0)

xvals = np.arange(0,A+1,1)

fit_params =  [[],[],[],[]]
fit = []
funcs = [fits.gaussian, fits.hyper_sec, fits.parabola, fits.lorentzian]
markers = ["x","o","v","s"]
colors = ["red","blue","green","black"]
labels = ["Gaussian", "Sech","Parabola", "Lorentzian"]

fig = plt.figure()

plt.bar(xvals, avg, align = "edge", width = 1,color="purple",edgecolor="black",zorder=0)
plt.xlabel("A")
plt.ylabel("Average Value of Casimir Operator")
plt.title("Amax="+str(A)+" Ntraj="+str(A**2))
z = 1

for i in range(len(fit_params)):
    fit_params[i], error = curve_fit(funcs[i],xvals,avg)
    fit.append(funcs[i](xvals,*fit_params[i]))
    plt.scatter(xvals,fit[i],color=colors[i],s=8,marker=markers[i],label=labels[i],zorder=z) 
    z += 1

plt.legend()
plt.savefig("/home/thomas/fluxtubes/figs/A"+str(A)+"fit.png")

