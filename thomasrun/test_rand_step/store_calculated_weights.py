import numpy as np
import os
import sys

Amax = int(sys.argv[1])

random_walks = np.loadtxt("/home/thomas/fluxtubes/thomasrun/test_rand_step/rand_traj/good_random_traj_"+str(Amax)+".txt", dtype = int, skiprows=1)
possible_a = np.array(range(1,Amax))
possible_p = np.array(range(np.amax(random_walks[:,1]+1)))
possible_q = np.array(range(np.amax(random_walks[:,2]+1)))
possible_n = np.array(range(np.amax(random_walks[:,3]+1)))
for a in possible_a:
    for p in possible_p:
        for q in possible_q:
            for n in possible_n:
                os.system("./calc_test_weight "+" "+str(a)+" "+str(p)+" "+str(q)+" "+str(n)+" "+str(Amax))
                random_walks = random_walks[(random_walks[:,0]!=a) & (random_walks[:,1]!=p) & (random_walks[:,2]!=q) & (random_walks[:,3]!=n),:]
