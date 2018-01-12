import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab

def massFlux_calc():
    """
    Inputs: none
    Return: none
    
    Takes in mdot data and prints out mass flux
    """
    
    massFluxNorm = 2.29434e-21
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Mdot','r')]
    
    xv1=[f1[i][0] for i in range(2,len(f1),288)]
    yv1=[f1[i][2] for i in range(2,len(f1),288)]
    
    mdot = 0
    for i in range(10,len(yv1)):
        print (yv1[i], mdot)
        mdot = mdot + yv1[i]*(xv1[i]-xv1[i-1])
    
    print (mdot/xv1[-1]/massFluxNorm)
    return 