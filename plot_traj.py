import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

whichA = sys.argv[1]
A = int(whichA)
filename = "A" +whichA+ ".txt"

traj_total = np.loadtxt("/home/thomas/fluxtubes/thomasrun/trajectories/" + filename)
average = np.average(traj_total, axis=0)
A_vals = np.arange(0,A+1,1)
x_vals = np.arange(0,A+2,2)


fig = plt.figure()
plt.bar(A_vals, average, align = "edge", width = 1,color="red",edgecolor="black")
plt.xlabel("A")
plt.xticks(x_vals)
plt.ylabel("Average Value of Casimir Operator")
plt.title("Amax="+str(A)+", Ntraj ="+str(A**2))
plt.savefig("/home/thomas/fluxtubes/figs/A"+whichA+".png")
