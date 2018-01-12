import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab

def MovingAverage(array, windowsize):
    average = []
    for i in range(len(array)):
        if i == len(array)/windowsize:
            break
        asum = 0
        lim2 = (i+1) * windowsize
        lim1 = i * windowsize
        for k in range(lim1,lim2):
            asum += array[k]
        asum /= windowsize
        average.append(asum)
    return average


def massFluxXmin_Tave():
    """
    Inputs: none
    Return: none
    
    Takes in mdot data and outputs figure labeled as 'massFluxXmin_Tave.eps'
    """
    massFluxNorm = 2.29434e-21

    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Mdot','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Mdot','r')]
    f3 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/Mdot','r')]
    
    xv1=[f1[i][0] for i in range(2,len(f1),256)]
    yv1=[f1[i][2] for i in range(2,len(f1),256)]
    
    xv2=[f2[i][0] for i in range(2,len(f2),256)]
    yv2=[f2[i][2] for i in range(2,len(f2),256)]
    
    xv3=[f3[i][0] for i in range(2,len(f3),256)]
    yv3=[f3[i][2] for i in range(2,len(f3),256)]
    
    xv1=np.asarray(xv1)/1e4
    xv2=np.asarray(xv2)/1e4
    xv3=np.asarray(xv3)/1e4
    
    yv1=abs(np.asarray(yv1))/massFluxNorm
    yv2=abs(np.asarray(yv2))/massFluxNorm
    yv3=abs(np.asarray(yv3))/massFluxNorm
    
    MAx1 = MovingAverage(xv1,10)
    MAx2 = MovingAverage(xv2,10)
    MAx3 = MovingAverage(xv3,10)
    MAy1 = MovingAverage(yv1,10)
    MAy2 = MovingAverage(yv2,10)
    MAy3 = MovingAverage(yv3,10)
    
    x1 = [0,9]
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'$\dot{m}$',**font)
    plt.xlabel(r'$t\,[\times 10^4 GM/c^3]$',**font)
    #pylab.xlim([0,1])
    pylab.ylim([1e-3,2e1])
    plt.yscale('log')
    plot1, = plt.plot(MAx1, MAy1, "k--", dashes=[6,4], color = 'red', mec = 'red', markersize = 6, linewidth = 2)
    plot2, = plt.plot(MAx2, MAy2, "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 2)
    plot3, = plt.plot(MAx3, MAy3, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    plot4, = plt.plot(x1,[0.01,0.01], "k--", dashes=[6,4], color = 'red', mec = 'red', markersize = 6, linewidth = 1)
    plot5, = plt.plot(x1,[1,1], "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 1)
    plot6, = plt.plot(x1,[10,10], "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 1)
    plt.legend((r'S01E',r'S1E',r'S10E'), 'upper right', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    #leg.get_frame().set_alpha(0)
    frame = leg.get_frame()
    frame.set_facecolor('white')
    frame.set_linewidth(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../figures/massFluxXmin_Tave.eps', format='eps', transparent='True')
    plt.show()
    return

