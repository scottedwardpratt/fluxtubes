import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os
import sys

Amax_list = [5,10,15,20]
max_dist = []
for Amax in Amax_list:
    os.system("./gluon_walker "+str(Amax))

    traj = np.loadtxt("gluon_walks/gluon_traj_"+str(Amax)+".txt", skiprows=1)
    
    for i in range(0,traj.shape[0],Amax+1):
        curr_traj = np.zeros((Amax +1,2))
        for j in range(Amax+1):
            curr_traj[j,0] = traj[i+j,1]
            curr_traj[j,1] = traj[i+j,2]
        dist = np.sqrt(curr_traj[:,0]**2 + curr_traj[:,1]**2)
        maxarg = np.argmax(dist)
        max_dist.append([curr_traj[maxarg,0],curr_traj[maxarg,1],Amax])

max_dist_array = np.array(max_dist)

fig = plt.figure(0)
ax = fig.add_subplot(111, aspect='equal')
ax.set_xlabel("p")
ax.set_ylabel("q")
ax.set_xlim(0,20)
ax.set_ylim(0,20)
ax.set_title("Furthest Points Reached by Returning Random Walks with Different Values of Amax")
for Amax in Amax_list:
    p = max_dist_array[(max_dist_array[:,2]==Amax),0]
    q = max_dist_array[(max_dist_array[:,2]==Amax),1]
    p_width = np.average(p)
    q_height = np.average(q)
    el = matplotlib.patches.Ellipse((0,0), p_width, q_height, angle=45, color=(np.random.random(),np.random.random(),np.random.random()), label="Amax="+str(Amax))
    ax.add_artist(el)
plt.show()
ax.legend()
plt.savefig("plots/contour.png")
