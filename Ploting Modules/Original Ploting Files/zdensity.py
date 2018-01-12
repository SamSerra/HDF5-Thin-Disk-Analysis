import sys
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def zdensity(file1, file2, file3, nt, nx, nskip, out):
    """
    Inputs:
        file1: filepath
        file2: filepath
        file3: filepath
        nt: int 
        nx: int
        nskip: int
        out: filepath
        
    Return: none
    
    Takes in data files, returns and saves plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[7] correspond to inputs file1, ..., out.
    """
    """
    if __name__ == "__main__":       #handels case where function called from command line
        tgas(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])
    """
    
    dunit = 1.409e+16
    f1 = open(file1,'r')
    f2 = open(file2,'r')
    f3 = open(file3,'r')
    
    x1f = []
    y1f = []
    z1f = []
    x2f = []
    y2f = []
    x3f = []
    y3f = []
    for line in f1:
        p1 = line.split()
        x1f.append(float(p1[0]))
        y1f.append(float(p1[1]))
        z1f.append(np.log10(np.abs(float(p1[2])*dunit)))
    for line in f2:
        p2 = line.split()
        x2f.append(float(p2[0]))
        y2f.append(float(p2[1]))
    for line in f3:
        p3 = line.split()
        x3f.append(float(p3[0]))
        y3f.append(float(p3[1]))
    
    x1 = []
    y1 = []
    z1 = []
    x2 = []
    y2 = []
    x3 = []
    y3 = []
    
    zmin = 100
    zmax = 0
    for j in range(nx):
        z1.append([])
        for i in range(int(nt/nskip)):
            if i == 0:
                y1.append(y1f[j])
            if j == 0:
                x1.append(x1f[nskip*i*nx]/1e4)
            z1[j].append(z1f[nskip*i*nx+j])
            if z1[j][i] < zmin:
                zmin = z1[j][i]
            if z1[j][i] > zmax:
                zmax = z1[j][i]
    
    levels = MaxNLocator(nbins=100).bin_boundaries(-10,-2)
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.ylabel(r'$\theta$',**font)                          #draws x-axis label
    plt.xlabel(r'$t\,[GM/c^3]$',**font)
    pylab.xlim([0,np.max(x1)])
    pylab.ylim([np.min(y1),np.max(y1)])
    #plt.xticks(np.arange(min(x1), max(x1)+1, 1000))
    # contours are *point* based plots, so convert our bound into point
    # centers
    plt.contourf(x1, y1, z1, levels=levels, extend="both")
    #plt.plot(x2f,y2f,'k',lw=2.5)
    #plt.plot(x3f,y3f,'k',lw=2.5)
    cbar = plt.colorbar(format="%.2f")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14)
    cbar.set_label(r'$\mathrm{log} \,\rho(z)\,[\mathrm{cgs}]$',**font)
    plt.savefig(out+".eps", format='eps', transparent='True')
    plt.show()
    return