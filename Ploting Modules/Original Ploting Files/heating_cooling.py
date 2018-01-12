import sys
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def heating_cooling(file1, file2, nt, nx, nskip, xmin, xmax, out):
    """
    Inputs: file1: filepath
            file 2: filepath
            nt: int
            nx: int
            nskip: int 
            xmin: float
            xmax: float
            out: str (file path and name, e.g. ../figures/genericname)
    Return: none
    
    Takes in 2 data files, returns and saves contour plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[8] correspond to inputs file1, ..., out.
    """
    """
    if __name__ == "__main__": #handles case where function called from command line
        file1 = sys.argv[1]
        file2 = sys.argv[2]
        nt = int(sys.argv[3])
        nx = int(sys.argv[4])
        nskip = int(sys.argv[5])
        xmin = float(sys.argv[6])
        xmax = float(sys.argv[7])
        out = sys.argv[8]
    """
    alpha = 0.02
    
    f1 = open(file1,'r')
    f2 = open(file2,'r')
    xf = []
    yf = []
    zf = []
    z1f = []
    
    for line in f1:
        p1 = line.split()
        xf.append(float(p1[1]))
        yf.append(float(p1[0]))
        zf.append(np.abs(float(p1[2])))
    for line in f2:
        p2 = line.split()
        z1f.append(np.abs(float(p2[2])))
    
    x = []
    y = []
    z1 = []
    tc = []
    
    
    zmin = 100
    zmax = 0
    for j in range(int(nt/nskip)):
        z1.append([])
        for i in range(nx):
            if i == 0:
                y.append(yf[nskip*j*nx]/1e4)
            if j == 0:
                x.append(xf[i])
                tc.append(6.28*(xf[i]**(1.5))/alpha/1e4)
            z1[j].append(np.log10(zf[nskip*j*nx+i]/z1f[nskip*j*nx+i]))
            if z1[j][i] < zmin:
                zmin = z1[j][i]
            if z1[j][i] > zmax:
                zmax = z1[j][i]
    
    level = MaxNLocator(nbins=100).bin_boundaries(-1.5,1.5)
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.xlabel(r'$R\,[GM/c^2]$',**font)                          #draws x-axis label
    plt.ylabel(r'$t\,[\times 10^4 GM/c^3]$',**font)             #draws y-axis label
    #plt.xscale('log')
    pylab.xlim([xmin,xmax])
    pylab.ylim([0,y[len(y)-1]])
    # contours are *point* based plots, so convert our bound into point
    # centers6
    plt.plot(x,tc,'w',lw=2)
    plt.savefig(out+".eps",format='eps',transparent='True')
    plt.show()
    
    plt.contourf(x, y, z1, levels=level, extend="both")
    cbar = plt.colorbar(format="%.2f")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14) 
    cbar.set_label(r'$\log\,(Q^+/Q^-)$',**font)
    plt.show()
    
    return 

