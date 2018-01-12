import sys
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def vertStress(file1, nt, nx, nskip, xmin,xmax,vmin,vmax,out):
    """
    Inputs:
        file1: filepath
        nt: int 
        nx: int
        nskip: int
        xmin: float
        xmax: float
        vmin: float
        vmax: float
        out: filepath
        
    Return: none
    
    Takes in data files, returns and saves plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[9] correspond to inputs file1, ..., out.
    """

    if __name__ == "__main__":       #handels case where function called from command line
        tgas(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7], sys.argv[8],sys.argv[9])

    alpha = 0.02
    f1 = open(file1,'r')
    
    
    xf = []     
    yf = []
    zf = []
    for line in f1:
        p1 = line.split()
        xf.append(float(p1[1]))
        yf.append(float(p1[0]))
        zf.append(float(p1[2]))
    
    x = []
    y = []
    z = []
    tLE = []
    
    zmin = 100
    zmax = -100
    for j in range(int(nt/nskip)):
        z.append([])
        for i in range(nx):
            if i == 0:
                y.append(yf[nskip*j*nx]/1.e4)
            if j == 0:
                x.append(xf[i])
                tLE.append(6.28*(xf[i]**(1.5))/alpha/1.e4)
            z[j].append(np.log10(zf[nskip*j*nx+i]))
            if z[j][i] < zmin:
                zmin = z[j][i]
            if z[j][i] > zmax:
                zmax = z[j][i]
    
    levels = MaxNLocator(nbins=100).bin_boundaries(vmin,vmax)
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.xlabel(r'$R\,[GM/c^2]$',**font)                          #draws x-axis label
    plt.ylabel(r'$t\,[\times 10^4 GM/c^3]$',**font)             #draws y-axis label
    pylab.xlim([xmin,xmax])
    pylab.ylim([0,y[len(y)-1]])
    # contours are *point* based plots, so convert our bound into point
    # centers
    plt.contourf(x, y, z, levels=levels, extend="both")
    plt.plot(x,tLE,'w',lw=2)
    cbar = plt.colorbar(format="%.1f")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14)
    cbar.set_label(r'$\log\,(W_{r\phi})\,[\mathrm{cgs}]$',**font)
    plt.savefig(out+".eps",format='eps',transparent='True')
    plt.show()
    return