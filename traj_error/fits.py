import numpy as np

def gaussian(x,a,b,c):
    return a*(np.exp(-b*(x-c)**2))

def hyper_sec(x,a,b,c):
    return a*(1/np.cosh(b*(x-c)))

def parabola(x,a,b,c):
    return a*x**2 + b*x +c

def lorentzian(x,a,b,c):
    return a*((0.5*b)/((x-c)**2 + (0.5*b)**2))

