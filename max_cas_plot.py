import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit


def linear_fit(x,m,b):
    return (x*m+b)

A = np.arange(10,101,1)
max_cas = []

for a in A:
    traj_total = np.loadtxt("/home/thomas/fluxtubes/thomasrun/trajectories/A" + str(a)+".txt")
    max_cas.append(np.amax(traj_total))

max_cas_array = np.array(max_cas)
params, error = curve_fit(linear_fit,A,max_cas_array)
fit = params[0]*A + params[1]


fig = plt.figure()
plt.xlabel("Amax")
plt.ylabel("Maximum value of Casimir")
plt.scatter(A,max_cas,color="red",label="Max Casimir")
plt.plot(A,fit,color="black",label="C = "+str(np.round(params[0],2))+"*A + "+ str(np.round(params[1],2)))
plt.title("Maximum Casimir Value vs. Amax")
plt.legend()
plt.savefig("/home/thomas/fluxtubes/figs/max_casimir.png")
