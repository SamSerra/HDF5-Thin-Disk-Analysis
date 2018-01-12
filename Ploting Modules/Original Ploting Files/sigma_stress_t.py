import sys
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab

def sigma_stress_t(file1, file2, nx, tstart, ymax1,ymin2, ymax2, out):
    """
    Inputs:
        file1: filepath
        file2: filepath
        nx: int
        tstart: int
        ymax1: float
        ymax2: float
        out: filepath
        
    Return: none
    
    Takes in data files, returns and saves plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[7] correspond to inputs file1, ..., out.
    """
    """
    if __name__ == "__main__":       #handels case where function called from command line
        sigma_tau_t(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])
    """
    
    istart=tstart*nx
    istop=(tstart+1)*nx
    
    f1 = [np.array(line.split()).astype('float') for line in open(file1,'r')]
    f2 = [np.array(line.split()).astype('float') for line in open(file2,'r')]
    
    xv1=[f1[i][1] for i in range(istart,istop)]
    yv1=[f1[i][2] for i in range(istart,istop)]
    
    xv2=[f2[i][1] for i in range(istart,istop)]
    yv2=[f2[i][2]/1e17 for i in range(istart,istop)]
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    fig, ax0 = plt.subplots()
    ax0.set_ylim([0,ymax1])
    ax0.set_xlim([xv1[0],xv1[len(xv1)-1]])
    #ax0.set_yscale('log')
    ax0.set_ylabel(r'$\Sigma$ [cgs]',**font)
    ax0.set_xlabel(r'$r\,[GM/c^2]$',**font)
    #ax0.set_xticklabels(**font)
    #ax0.set_yticklabels(**font)
    ax0.tick_params(length=5, which='minor')
    ax0.tick_params(length=8, which='major')
    ax0.plot(xv1, yv1, "k--", dashes=[6,4], markersize = 6, linewidth = 2)
    #plot0.set_dashes([6,4])
    ax1 = ax0.twinx()
    ax1.set_ylabel(r'$W_{r\phi}$ [$\times 10^{17}$ cgs]',**font)
    ax1.set_ylim([ymin2,ymax2])
    ax1.set_xlim([xv2[0],25])
    #ax1.set_yscale('log')
    ax1.tick_params(length=5, which='minor')
    ax1.tick_params(length=8, which='major')
    axtext = ax1.yaxis.get_offset_text()
    plt.setp(axtext,**font)
    ax1.plot(xv2, yv2, "k--", dashes=[2,3], color = 'blue', mec = 'blue', markersize = 6, linewidth = 2)
    
    
    #plt.legend((r'$\Sigma$',r'$W_{r,\phi}$'), 'lower right', shadow=False, handlelength=3)
    #leg =  plt.gca().get_legend()
    #leg.get_frame().set_alpha(0)
    #frame = leg.get_frame()
    #frame.set_facecolor('white')
    #frame.set_linewidth(0)
    #ltext = leg.get_texts()
    #plt.setp(ltext,**font)
    axtext = ax0.get_xmajorticklabels()
    plt.setp(axtext,**font)
    axtext = ax0.get_ymajorticklabels()
    plt.setp(axtext,**font)
    axtext = ax1.get_ymajorticklabels()
    plt.setp(axtext,**font)
    
    plt.savefig(out+".eps", format='eps', transparent='True')
    plt.show()
    return