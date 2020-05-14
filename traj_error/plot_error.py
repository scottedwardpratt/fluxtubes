import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

Amax = int(sys.argv[1])
log = sys.argv[2]
error = np.loadtxt("/home/thomas/fluxtubes/traj_error/error_data_"+str(Amax)+".txt")

x_ticks = (10, Amax, 10)


if log == "log":
    fig = plt.figure()

    plt.xlabel("Amax")
    plt.ylabel("Log10 Mean Squared Error")

    plt.scatter(error[:,0],np.log10(error[:,1]),color="red",marker="x",label="Gaussian")
    plt.scatter(error[:,0],np.log10(error[:,2]),color="blue",marker="o",label="Sech")
    plt.scatter(error[:,0],np.log10(error[:,3]),color="green",marker="v",label="Parabola")
    plt.scatter(error[:,0],np.log10(error[:,4]),color="black",marker="s",label="Lorentzian")

    plt.title("Lob10 MSE vs. Amax for Different Fits")

    plt.legend()

    plt.savefig("/home/thomas/fluxtubes/figs/traj_fit_log_error.png")

elif log == "notlog":
    fig = plt.figure()

    plt.xlabel("Amax")
    plt.ylabel("Mean Squared Error")

    plt.scatter(error[:,0],error[:,1],color="red",marker="x",label="Gaussian")
    plt.scatter(error[:,0],error[:,2],color="blue",marker="o",label="Sech")
    plt.scatter(error[:,0],error[:,3],color="green",marker="v",label="Parabola")
    plt.scatter(error[:,0],error[:,4],color="black",marker="s",label="Lorentzian")

    plt.title("MSE vs. Amax for Different Fits")

    plt.legend()

    plt.savefig("/home/thomas/fluxtubes/figs/traj_fit_error.png")
else:
    print("Please enter log or notlog as an argument")

