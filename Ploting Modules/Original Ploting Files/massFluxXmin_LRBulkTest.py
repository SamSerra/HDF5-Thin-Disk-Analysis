import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab


def massFluxXmin_LRBulkTest():
    """
    Input: none
    Return: none
    
    Takes in mdot data and returns figure labeled as 'massFluxXmin.eps'
    """
    massFluxNorm = 2.29434e-21
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_lr_PP/Mdot','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_NoVisc_PP/Mdot','r')]
    f3 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_bulk_PP/Mdot','r')]
    f4 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_nppm_PP/Mdot','r')]
    f5 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_zmin0_PP/Mdot','r')]
    
    xv1=[f1[i][0] for i in range(2,len(f1),144)]
    yv1=[f1[i][2] for i in range(2,len(f1),144)]
    
    xv2=[f2[i][0] for i in range(2,len(f2),144)]
    yv2=[f2[i][2] for i in range(2,len(f2),144)]
    
    xv3=[f3[i][0] for i in range(2,len(f3),144)]
    yv3=[f3[i][2] for i in range(2,len(f3),144)]
    
    xv4=[f4[i][0] for i in range(2,len(f4),144)]
    yv4=[f4[i][2] for i in range(2,len(f4),144)]
    
    xv5=[f5[i][0] for i in range(2,len(f5),144)]
    yv5=[f5[i][2] for i in range(2,len(f5),144)]
    
    xv1=np.asarray(xv1)/1e4
    xv2=np.asarray(xv2)/1e4
    xv3=np.asarray(xv3)/1e4
    xv4=np.asarray(xv4)/1e4
    xv5=np.asarray(xv5)/1e4
    
    yv1=abs(np.asarray(yv1))/massFluxNorm
    yv2=abs(np.asarray(yv2))/massFluxNorm
    yv3=abs(np.asarray(yv3))/massFluxNorm
    yv4=abs(np.asarray(yv4))/massFluxNorm
    yv5=abs(np.asarray(yv5))/massFluxNorm
    
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'$\dot{m}$',**font)
    plt.xlabel(r'$t\,[\times 10^4 GM/c^3]$',**font)
    pylab.xlim([0,1])
    pylab.ylim([2e-5,1000])
    plt.yscale('log')
    plot1, = plt.plot(xv1, yv1, "k--", dashes=[6,4], color = 'red', mec = 'red', markersize = 6, linewidth = 2)
    plot2, = plt.plot(xv2, yv2, "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 2)
    plot3, = plt.plot(xv3, yv3, "k--", dashes=[6,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    plot4, = plt.plot(xv4, yv4, "k--", dashes=[8,6,2,4,2,4], color = 'cyan', mec = 'cyan', markersize = 6, linewidth = 2)
    plot5, = plt.plot(xv5, yv5, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'blue', mec = 'blue', markersize = 6, linewidth = 2)
    plt.legend((r'S01E_lr',r'S01E_NoVisc',r'S01E_bulk',r'S01E_nppm',r'S01E_zmin0'), 'upper right', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    #plt.savefig('../figures/massFluxXmin.eps', format='eps', transparent='True')
    plt.show()
    return 
