import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

plt.figure(figsize = (10,10))
plt.xlabel("p")
plt.ylabel("q")
plt.title("Fartherst Point Reached by Returning Random Walks with Different Values of Amax")
plt.xlim(0,20)
plt.ylim(0,20)

for a in Amax_list:
    if a == 5:
        color = "red"
    if a ==10:
        color = "blue"
    if a ==15:
        color = "green"
    if a ==20:
        color = "black"
    plt.scatter(max_dist_array[(max_dist_array[:,2]==a),0], max_dist_array[(max_dist_array[:,2]==a),1], color=color,label="Amax="+str(a))

plt.legend()
plt.show()
plt.savefig("plots/contour.png")
