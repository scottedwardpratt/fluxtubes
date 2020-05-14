import numpy as np
import os
import sys
from scipy.optimize import curve_fit
import fits

def get_average(curr_A):
    traj_total = np.loadtxt("/home/thomas/fluxtubes/thomasrun/trajectories/A" + curr_A +".txt")
    return np.average(traj_total, axis=0)

def mean_squared_error(data,fit):
    return np.sum((fit-data)**2)/fit.size

Amax = int(sys.argv[1])
A_vals = np.arange(101,Amax+1,1)

os.chdir("/home/thomas/fluxtubes/thomasrun/")

param_error = [[],[],[],[]]
fit_params = [[],[],[],[]]
funcs = [fits.gaussian, fits.hyper_sec, fits.parabola, fits.lorentzian]

matrix = []

for A in A_vals:
    
    row = []
    fit = []
    fit_error = []
    
    row.append(A)

    os.system("./traj " + str(A))

    avg = get_average(str(A))
    xvals = np.arange(0,A+1,1)
    
    for i in range(len(fit_params)):
        fit_params[i], param_error[i] = curve_fit(funcs[i],xvals,avg)
        fit.append(funcs[i](xvals,*fit_params[i]))
        fit_error.append(mean_squared_error(avg,fit[i]))
        row.append(fit_error[i])
    
    matrix.append(row)

os.chdir("/home/thomas/fluxtubes/")

error_data = np.array(matrix)

np.savetxt("/home/thomas/fluxtubes/traj_error/error_data_"+str(Amax)+"_.txt",error_data)
