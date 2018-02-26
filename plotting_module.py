"""
Created: Fri Oct 20 16:30:05 2017
Author: samps

Bhupendra's Useful Module for Plotting Stuff (BUMPS)
Call functions with __name__.<function_name>()

"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import colors, ticker, cm
import scipy.fftpack
from scipy import signal
import pylab
import pylab, math
from matplotlib import rcParams
from math import sqrt
from math import cos
from math import acos
from math import pi
from math import log

rcParams.update({'figure.autolayout': True})

#%% Suplementary Functions

def MovingAverage(array, windowsize):
    """
    Inputs: array, windowsize
    Return: average6
    
    Takes in array and windowsize and returns array with moving average
    """
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

#%% efficiency.py

def efficiency():
    """
    Input: none
    Return: none
    
    Takes in Mdot and luminosity data, manipulates and plots them, and saves the figure as efficiency.eps
    """
    luminosityNorm = 2.29434e-21
    massFluxNorm = 2.29434e-21
    timeUnit = 3.26111e-05
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/luminosity','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Mdot','r')]
    f3 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/luminosity','r')]
    f4 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Mdot','r')]
    f5 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/luminosity','r')]
    f6 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/Mdot','r')]
    
    xv1=[f1[i][0] for i in range(0,len(f1))]
    yv1=[f1[i][1] for i in range(0,len(f1))]
        
    xv2=[f2[i][0] for i in range(2,len(f2),256)]
    yv2=[f2[i][2] for i in range(2,len(f2),256)]
        
    xv3=[f3[i][0] for i in range(0,len(f3))]
    yv3=[f3[i][1] for i in range(0,len(f3))]
    
    xv4=[f4[i][0] for i in range(2,len(f4),256)]
    yv4=[f4[i][2] for i in range(2,len(f4),256)]
        
    xv5=[f5[i][0] for i in range(0,len(f5))]
    yv5=[f5[i][1] for i in range(0,len(f5))]
        
    xv6=[f6[i][0] for i in range(2,len(f6),256)]
    yv6=[f6[i][2] for i in range(2,len(f6),256)]
    
    xv1=np.asarray(xv1)/1e4
    xv3=np.asarray(xv3)/1e4
    xv5=np.asarray(xv5)/1e4
    
    y1 = []
    y3 = []
    y5 = []
    for i in range(0,len(xv1)):
        y1.append(abs(yv1[i])/abs(yv2[i]))
    for i in range(0,len(xv3)):
        y3.append(abs(yv3[i])/abs(yv4[i]))
    for i in range(0,len(xv5)):
        y5.append(abs(yv5[i])/abs(yv6[i]))
    
    eta0 = [0.057 for i in range(0,len(f3))]
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'$\eta$',**font)
    plt.xlabel(r'$t\,[\times 10^4 GM/c^3]$',**font)
    #pylab.xlim([0,1])
    pylab.ylim([1e-4,10])
    plt.yscale('log')
    plot1, = plt.plot(xv1, y1, "k--", dashes=[6,4], color = 'red', mec = 'red', markersize = 6, linewidth = 2)
    plot2, = plt.plot(xv3, y3, "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 2)
    plot3, = plt.plot(xv5, y5, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    plot4, = plt.plot(xv3,eta0, 'k--', color = 'gray', mec = 'gray', markersize = 6, linewidth = 1)
    plt.legend((r'S01E',r'S1E',r'S10E'), 'upper right', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../figures/efficiency.eps', format='eps', transparent='True')
    plt.show()
    return

#%% efficiency_Tave.py

def efficiency_Tave(MovingAverage):
    """
    Inputs: MovingAverage
    Return: none
    
    Takes in Mdot and luminosity data, manipulates and plots it, and saves figure as 'efficiency_Tave.eps'
    """
    luminosityNorm = 2.29434e-21
    massFluxNorm = 2.29434e-21
    timeUnit = 3.26111e-05
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/luminosity','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Mdot','r')]
    f3 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/luminosity','r')]
    f4 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Mdot','r')]
    f5 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/luminosity','r')]
    f6 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/Mdot','r')]
    
    xv1=[f1[i][0] for i in range(0,len(f1))]
    yv1=[f1[i][1] for i in range(0,len(f1))]
        
    xv2=[f2[i][0] for i in range(2,len(f2),256)]
    yv2=[f2[i][2] for i in range(2,len(f2),256)]
        
    xv3=[f3[i][0] for i in range(0,len(f3))]
    yv3=[f3[i][1] for i in range(0,len(f3))]
    
    xv4=[f4[i][0] for i in range(2,len(f4),256)]
    yv4=[f4[i][2] for i in range(2,len(f4),256)]
        
    xv5=[f5[i][0] for i in range(0,len(f5))]
    yv5=[f5[i][1] for i in range(0,len(f5))]
        
    xv6=[f6[i][0] for i in range(2,len(f6),256)]
    yv6=[f6[i][2] for i in range(2,len(f6),256)]
    
    xv1=np.asarray(xv1)/1e4
    xv3=np.asarray(xv3)/1e4
    xv5=np.asarray(xv5)/1e4
    
    y1 = []
    y3 = []
    y5 = []
    for i in range(0,len(xv1)):
        y1.append(abs(yv1[i])/abs(yv2[i]))
    for i in range(0,len(xv3)):
        y3.append(abs(yv3[i])/abs(yv4[i]))
    for i in range(0,len(xv5)):
        y5.append(abs(yv5[i])/abs(yv6[i]))
    
    eta0 = [0.057 for i in range(0,len(f3))]
    
    MA1 = MovingAverage(y1,10)
    MA2 = MovingAverage(y3,10)
    MA3 = MovingAverage(y5,10)
    MAx1 = MovingAverage(xv1,10)
    MAx2 = MovingAverage(xv3,10)
    MAx3 = MovingAverage(xv5,10)
    
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'$\eta$',**font)
    plt.xlabel(r'$t\,[\times 10^4 GM/c^3]$',**font)
    #pylab.xlim([0,1])
    pylab.ylim([1e-3,10])
    plt.yscale('log')
    plot1, = plt.plot(MAx1, MA1, "k--", dashes=[6,4], color = 'red', mec = 'red', markersize = 6, linewidth = 2)
    plot2, = plt.plot(MAx2, MA2, "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 2)
    plot3, = plt.plot(MAx3, MA3, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    plot4, = plt.plot(xv3,eta0, 'k--', color = 'gray', mec = 'gray', markersize = 6, linewidth = 1)
    plt.legend((r'S01E',r'S1E',r'S10E'), 'upper right', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    #plt.plot(xv1,yv1)
    plt.savefig('../figures/efficiency_Tave.eps', format='eps', transparent='True')
    plt.show()
    return

#%% freq_rP.py

def freq_rP():
    """
    Inputs: none
    Return: none
    
    Imports Mdot data, preforms DFT, and plots figure
    """
    
    G     = 6.67259e-8
    Mbh   = 6.62*1.988e+33
    c     = 2.99792458e+10
    tunit = G*Mbh/c**3.
    font = {'fontname':'Times New Roman','fontsize':18, 'weight':'normal'}
    vmin = 5
    vmax = 12.
    f1 = [np.array(line.split()).astype('float')
    for line in open('../shakura_01E_128x96_PP/Mdot','r')]
    
    trange = len(f1)
    nr = 128
    xv1  = [f1[i][0] for i in range(0,trange,nr)]
    rv1  = [f1[i][1] for i in range(0,nr,1)]
    fftfreq1 = np.fft.rfftfreq(len(xv1), (xv1[1] - xv1[0]))
    
    fft1 = []
    r    = []
    for j in range(0,nr,5):
    	yv1  = [f1[i][2] for i in range(j,trange,nr)]
    	A    = np.fft.rfft(yv1)
    	fft1.append(2*np.abs(A)**2*(xv1[len(xv1)-1]-xv1[0])/np.abs(A[0])**2)
    	r.append(rv1[j])
    
    
    # fft1 = np.log10(fft1)
    fft1 = np.transpose(fft1)
    fftfreq1 = fftfreq1*1.e3
    levels = MaxNLocator(nbins=100).bin_boundaries(vmin,vmax)
    
    plt.xticks(**font)
    plt.yticks(**font)             #draws y-axis label
    
    plt.contourf(r,fftfreq1,fft1, levels=levels, extend="both", cmap=cm.hot)
    cbar = plt.colorbar(format="%.1e")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(16)
    cbar.set_label(r'$\mathrm{PDS}\,P_\mathrm{tot}$',**font)
    # cbar.set_label(r'$\log\,\mathrm{PDS}\,\dot{m}$',**font)
    rr  = np.linspace(6.,40.,500000)
    omegaz = np.sqrt(1./rr**3.)/2./np.pi*1.e3
    plt.tick_params(axis='x', direction='out')
    plt.tick_params(axis='y', direction='out')
    omegar = np.sqrt(1 - 6./rr)*omegaz #Schwarschild
    plt.plot(rr,omegar,'b',lw=1.5)
    plt.plot(rr,omegaz,'w--',lw=1.5)
    # plt.yscale('log')
    pylab.xlim([5.5,15])
    pylab.ylim([0,5])
    plt.xlabel(r'$R\,[GM/c^2]$',**font)                          #draws x-axis label
    plt.ylabel(r'$\nu\,[\times 10^{-3}c^3/GM]$',**font)
    #plt.savefig('../figures/PDSPtot0_01E.eps',format='eps',Transparent='True')
    
    # r1 = 6.5
    # r2 = 6.55
    # for j in range(0,256,1):
    # 	if f1[j][1] > r1 and f1[j][1] < r2 : rj = j
    
    # yv1  = [f1[i][2] for i in range(rj,len(f1),256)]
    # yv1 = np.asarray(yv1)*1.e4/yv1[0]
    # A    = np.fft.rfft(yv1)
    # fft1 = 2*np.abs(A)**2*(xv1[len(xv1)-1]-xv1[0])/np.abs(A[0])**2
    # plt.plot(fftfreq1,fft1)
    # pylab.xlim([0.002,0.01])
    # # pylab.ylim([0,30000])
    # plt.xlabel(u'$\\nu\,[c^3/GM]$',**font)                          #draws x-axis label
    # plt.ylabel(ur'$\mathrm{Power}$',**font)
    # plt.yscale('log')
    # plt.savefig('../figures/PDSmdot01E_r8.eps',format='eps',Transparent='True')
    plt.show()


#%% gmode.py
    
def gmode():
    """
    Input: none
    Return: none
    
    Takes in Mdot data, plots PDSMdot figure, and saves it as 'PDSmdot01E_r7p7_8p0_8p3.eps'
    """
    G     = 6.67259e-8
    Mbh   = 6.62*1.988e+33
    c     = 2.99792458e+10
    tunit = G*Mbh/c**3.
    
    # f1 = [np.array(line.split()).astype('float') 
    # for line in open('../shakura_01E_PP/Ptot0','r')]
    # f2 = [np.array(line.split()).astype('float') 
    # for line in open('../shakura_1E_PP/Ptot0','r')]
    # f3 = [np.array(line.split()).astype('float') 
    # for line in open('../shakura_10E_PP/Ptot0','r')]
    
    f1 = [np.array(line.split()).astype('float') 
    for line in open('../shakura_01E_PP/Mdot','r')]
    f2 = [np.array(line.split()).astype('float')  
    for line in open('../shakura_1E_PP/Mdot','r')]
    f3 = [np.array(line.split()).astype('float')  
    for line in open('../shakura_10E_PP/Mdot','r')]
    
    
    r1 = 7.9
    r2 = 8.01
    
    r11 = 7.65
    r21 = 7.78
    
    r12 = 8.28
    r22 = 8.35
    
    r13 = 8.9
    r23 = 9.0
    
    
    for j in range(0,256,1):
    	if f1[j][1] > r1 and f1[j][1] < r2 : rj1 = j
    	if f2[j][1] > r1 and f2[j][1] < r2 : rj2 = j
    	if f3[j][1] > r1 and f3[j][1] < r2 : rj3 = j
    	
    	if f1[j][1] > r11 and f1[j][1] < r21 : rj11 = j
    	if f1[j][1] > r12 and f1[j][1] < r22 : rj12 = j
    	if f1[j][1] > r13 and f1[j][1] < r23 : rj13 = j
    
    	if f2[j][1] > r11 and f2[j][1] < r21 : rj21 = j
    
    	if f3[j][1] > r11 and f3[j][1] < r21 : rj31 = j
    	if f3[j][1] > r12 and f3[j][1] < r22 : rj22 = j
    	if f3[j][1] > r13 and f3[j][1] < r23 : rj23 = j
    
    print (f1[rj11][1],f1[rj12][1])
    
    xv1  = [f1[i][0] for i in range(0,len(f1),256)]
    yv1  = [f1[i][2] for i in range(rj1,len(f1),256)]
    
    yv11 = [f1[i][2] for i in range(rj11,len(f1),256)]
    yv12 = [f1[i][2] for i in range(rj12,len(f1),256)]
    yv13 = [f1[i][2] for i in range(rj13,len(f1),256)]
    
    xv2  = [f2[i][0] for i in range(0,len(f2),256)]
    yv2  = [f2[i][2] for i in range(rj2,len(f2),256)]
    
    xv21 = [f2[i][0] for i in range(0,len(f2),256)]
    yv21 = [f2[i][2] for i in range(rj21,len(f2),256)]
    yv22 = [f2[i][2] for i in range(rj22,len(f2),256)]
    yv23 = [f2[i][2] for i in range(rj23,len(f2),256)]
    
    xv3  = [f3[i][0] for i in range(0,len(f3),256)]
    yv3  = [f3[i][2] for i in range(rj3,len(f3),256)]
    
    xv31 = [f3[i][0] for i in range(0,len(f3),256)]
    yv31 = [f3[i][2] for i in range(rj31,len(f3),256)]
    
    
    xv1 = np.asarray(xv1)
    
    # yv1 = abs(np.asarray(yv1)*1.e4/yv1[0])
    # yv11= abs(np.asarray(yv11)*1.e4/yv11[0])
    # yv12= abs(np.asarray(yv12)*1.e4/yv12[0])
    # yv13= abs(np.asarray(yv13)*1.e4/yv13[0])
    # dt1 = (xv1[1] - xv1[0])
    # n1  = len(xv1)
    
    # xv2 = np.asarray(xv2)
    # xv21= np.asarray(xv21)
    # yv2 = abs(np.asarray(yv2)*1.e4/yv2[0])
    # yv21= abs(np.asarray(yv21)*1.e4/yv21[0])
    # yv22= abs(np.asarray(yv22)*1.e4/yv22[0])
    # yv23= abs(np.asarray(yv23)*1.e4/yv23[0])
    # dt2 = (xv2[1] - xv2[0])
    # n2  = len(xv2)
    
    # xv3 = np.asarray(xv3)
    # xv31= np.asarray(xv31)
    # yv3 = abs(np.asarray(yv3)*1.e4/yv3[0])
    # yv31= abs(np.asarray(yv31)*1.e4/yv31[0])
    # dt3 = (xv3[1] - xv3[0])
    # n3  = len(xv3) 
    
    # plt.psd(yv1,1000,1/dt1,scale_by_freq=None,c='r',label='7',ls='solid',lw=2)
    # plt.psd(yv11,1000,1/dt1,scale_by_freq=None,c='g',label='7.7',ls='dashed',lw=2)
    # plt.psd(yv12,1000,1/dt1,scale_by_freq=None,c='b',label='8.3',ls='dotted',lw=2)
    # plt.psd(yv13,1000,1/dt1,scale_by_freq=None,c='r',label='12.0',ls='solid',lw=2)
    
    # plt.psd(yv2,n2,1/dt2,scale_by_freq=None,c='g',label='6.0',ls='dashed',lw=1.5)
    # plt.psd(yv21,n2,1/dt2,scale_by_freq=None,c='g',label='7.0',ls='dashed',lw=2)
    # plt.psd(yv2,n2,1/dt2,scale_by_freq=None,c='g',label='1E',ls='solid',lw=2)
    # plt.psd(yv23,n2,1/dt2,scale_by_freq=None,c='b',label='9.0',ls='solid',lw=1.5)
    
    # plt.psd(yv3,n3,1/dt3,scale_by_freq=None,ls='solid', c = 'orange',label='10E',lw=2)
    # plt.psd(yv31,n3,1/dt3,scale_by_freq=None,ls='solid', c = 'orange',label='8.0',lw=2)
    
    # plt.text(0.0036,50, r'$\nu_\mathrm{g}$', fontsize=22)
    # plt.grid(False)
    # plt.gca().ticklabel_format(style='sci', axis='x')
    # font = {'fontname':'Times New Roman','fontsize':18, 'weight':'normal'}
    # plt.xticks(**font)
    # plt.yticks(**font)
    # plt.minorticks_on()
    # # plt.xscale('log')
    # # plt.yscale('log')
    # plt.xlim([0.002,0.016])
    # plt.ylim([46,140])
    plt.tick_params(length=5,which='minor')
    plt.tick_params(length=10,which='major')
    # plt.xlabel(u'$\\nu\,[c^3/GM]$',**font)                     
    # plt.ylabel(u'$\mathrm{PDS}\,[\dot{m}]$',**font)
    # plt.legend(loc=1,ncol=4,prop={'size':20},frameon=False)
    # plt.savefig('../figures/PDSmdot01E_r7p7_8p3.eps',format='eps',Transparent='True')
    
    font = {'fontname':'Times New Roman','fontsize':18, 'weight':'normal'}
    fftfreq1 = np.fft.rfftfreq(len(xv1), (xv1[1] - xv1[0]))
    fftfreq1 = fftfreq1*1.e3
    A1    = np.fft.rfft(yv11/yv11[0])
    A2    = np.fft.rfft(yv1/yv1[0])
    A3    = np.fft.rfft(yv12/yv12[0])
    fft1 = 2*np.abs(A1)**2*(xv1[len(xv1)-1]-xv1[0])/np.abs(A1[0])**2
    fft2 = 2*np.abs(A2)**2*(xv1[len(xv1)-1]-xv1[0])/np.abs(A2[0])**2
    fft3 = 2*np.abs(A3)**2*(xv1[len(xv1)-1]-xv1[0])/np.abs(A3[0])**2
    plt.plot(fftfreq1,fft1,color='g',label='7.7',ls='dashed',lw=2)
    plt.plot(fftfreq1,fft2,color='r',label='8.0',ls='dotted',lw=2)
    plt.plot(fftfreq1,fft3,color='b',label='8.3',ls='solid',lw=2)
    
    plt.legend(loc=1,ncol=4,prop={'size':20},frameon=False)
    plt.vlines(3.5,0.01,1e8,color='k',linestyle='dotted',lw=1)
    plt.vlines(10.67,0.1,1e5,color='b',linestyle='dashed',lw=1)
    plt.vlines(11.37,0.1,1e5,color='r',linestyle='dashed',lw=1)
    plt.vlines(11.8,0.1,1e5,color='g',linestyle='dashed',lw=1)
    
    plt.text(3,10, r'$3.50$', fontsize=16,color='k',rotation=90)
    plt.text(0.7,2e4, r'$\mathrm{01E}$', fontsize=18,color='k',rotation=0)
    # plt.text(10.2,38, r'$10.67$', fontsize=16,color='b',rotation=90)
    # plt.text(11,38, r'$11.37$', fontsize=16,color='r',rotation=90)
    # plt.text(11.9,38, r'$11.80$', fontsize=16,color='g',rotation=90)
    
    pylab.xlim([0,10])
    pylab.ylim([1e-1,1.e5])
    plt.xlabel(r'$\\nu\,[\\times 10^{-3}c^3/GM]$',**font)                          #draws x-axis label
    plt.ylabel(r'$\mathrm{PDS}\,\dot{m}$',**font)
    # plt.ylabel(ur'$\mathrm{PDS}\,P_\mathrm{tot}$',**font)
    plt.yscale('log')
    plt.xticks(**font)
    plt.yticks(**font) 
    
    a = plt.axes([.65, .65, .25, .2])
    a = plt.plot(fftfreq1,fft1,color='g',label='7.7',ls='dashed',lw=2)
    a = plt.plot(fftfreq1,fft2,color='r',label='8.0',ls='dotted',lw=2)
    a = plt.plot(fftfreq1,fft3,color='b',label='8.3',ls='solid',lw=2)
    a = plt.xlim([3.45,3.65])
    a = plt.ylim([11000,35000])
    # a = plt.yscale('log')
    a = plt.tick_params(length=5,which='minor')
    a = plt.tick_params(length=10,which='major')
    
    plt.savefig('../figures/PDSmdot01E_r7p7_8p0_8p3.eps',format='eps',Transparent='True')
    plt.show()
    return

#%% gmode_cavity.py
    
def gmode_cavity ():
    """
    Input: none
    Return: none
    
    Plots gmode cavity and saves figure as .eps
    """
    # G = 6.67e-8
    # c = 2.9999e+10
    # msun = 1.989e+33
    # mbh  = 6.62
    # munit = mbh
    # lunit = G*mbh/c**2
    # tunit = lunit/c
    # dunit = munit/lunit/lunit/lunit
    # edunit = 1.4e+37
    # medd   = 1.39e+17*6.62 #cgs units
    
    rr  = np.linspace(6.,20.,500)
    #omegar = np.sqrt((rr - 6)/rr/(rr - 2)**3)/2./np.pi #PW
    omegaz = np.sqrt(1./rr**3.)
    omegar = np.sqrt(1 - 6./rr)*omegaz/2./np.pi*1.e3 #Schwarschild
    
    plt.plot(rr,omegar,'b',lw=1.5)
    plt.hlines(3.5, 7.7,8.34,linewidth=3, color='k',linestyle='dotted')
    plt.vlines(12.2, 0.,3.5, linewidth=2, color='g',linestyle='dashed')
    plt.vlines(6.25, 0.,3.5,linewidth=2, color='g',linestyle='dashed')
    plt.text(7.1,5.6e-3,'$\\kappa_\mathrm{max} = 5.517\\times 10^{-3}$',fontsize=18)
    plt.text(7.1,3.53e-3,'$\\kappa = 3.507\\times 10^{-3}$',fontsize=18)
    font = {'fontname':'Times New Roman','fontsize':20, 'weight':'normal'}
    
    plt.xticks(**font)
    plt.yticks(**font)
    plt.xlabel(r'$r\,[GM/c^2]$',**font)                          #draws x-axis label
    plt.ylabel(r'$\\kappa\,[\\times 10^{-3}\,c^3/GM]$',**font)             #draws y-axis label
    #plt.xscale('log')
    pylab.xlim([7.6,8.4])
    pylab.ylim([3.47,3.52]) #RADP
    plt.savefig('../figures/gmode_cavity.eps',format='eps',Transparent='True')
    plt.show()
    return 

#%% heating_cooling.py
    
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
    if __name__ == "__main__": #handles case where function called from command line
        heating_cooling(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8])

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
    # centers
    plt.contourf(x, y, z1, levels=level, extend="both")
    plt.plot(x,tc,'w',lw=2)          
    cbar = plt.colorbar(format="%.2f")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14)
    cbar.set_label(r'$\log\,(Q^+/Q^-)$',**font)
    plt.savefig(out+".eps",format='eps',transparent='True')
    plt.show()
    return 
    
#%% height.py
    
def height(file1, nt, nx, nskip, xmin, xmax, vmin, vmax, out):
    """
    Inputs: file1: filepath
            nt: int
            nx: int
            nskip: int 
            xmin: float
            xmax: float
            vmin: float
            vmax: float
            out: str (file path and name, e.g. ../figures/genericname)
    Return: none
    
    Takes in data files, returns and saves contour plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[9] correspond to inputs file1, ..., out.
    """
    if __name__ == "__main__": #handles case where function called from command line
        height(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9])
    
    f1 = open(file1,'r')
    alpha = 0.02
    
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
    tc = []
    
    zmin = 100
    zmax = 0
    for j in range(int(nt/nskip)):
        z.append([])
        for i in range(nx):
            if i == 0:
                y.append(yf[nskip*j*nx]/1e4)
            if j == 0:
                x.append(xf[i])
                tc.append(6.28*(xf[i]**(1.5))/alpha/1e4)
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
    plt.plot(x,tc,'w',lw=2)
    cbar = plt.colorbar(format="%.3f")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14)
    # cbar.set_label(r'$H/R$',**font)
    # cbar.set_label(r'$H\,[GM/c^2]$',**font)
    cbar.set_label(r'$\log\,[\langle H \rangle/(GM/c^2)]$',**font)
    plt.savefig(out+".eps",format='eps',transparent='True')
    plt.show()
    return

#%% luminosity_Tave.py
    
def luminosity_Tave():
    """
    Inputs: none
    Return: none
    
    Takes in luminosity data and returns many plots
    """
    
    luminosityNorm = 2.29434e-21
    timeUnit = 3.26111e-05
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP_new/luminosity','r')]
    f1a = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP_new/luminosity_clean','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/luminosity','r')]
    f3 = [np.array(line.split()).astype('float') for line in open('../shakura_3E_PP/luminosity','r')]
    f4 = [np.array(line.split()).astype('float') for line in open('../shakura_3Ep_PP/luminosity','r')]
    f5 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP_new/luminosity','r')]
    
    xv1=[f1[i][0] for i in range(0,len(f1))]
    yv1=[f1[i][1] for i in range(0,len(f1))]
    
    xv1a=[f1a[i][0] for i in range(0,len(f1a))]
    yv1a=[f1a[i][1] for i in range(0,len(f1a))]
    
    xv2=[f2[i][0] for i in range(0,len(f2))]
    yv2=[f2[i][1] for i in range(0,len(f2))]
    
    xv3=[f3[i][0] for i in range(0,len(f3))]
    yv3=[f3[i][1] for i in range(0,len(f3))]
    
    xv4=[f4[i][0] for i in range(0,len(f4))]
    yv4=[f4[i][1] for i in range(0,len(f4))]
    
    xv5=[f5[i][0] for i in range(0,len(f5))]
    yv5=[f5[i][1] for i in range(0,len(f5))]
    
    xv1=np.asarray(xv1)/1e4
    xv1a=np.asarray(xv1a)/1e4
    xv2=np.asarray(xv2)/1e4
    xv3=np.asarray(xv3)/1e4
    xv4=np.asarray(xv4)/1e4
    xv5=np.asarray(xv5)/1e4
    
    yv1=abs(np.asarray(yv1))/luminosityNorm
    yv1a=abs(np.asarray(yv1a))/luminosityNorm
    yv2=abs(np.asarray(yv2))/luminosityNorm
    yv3=abs(np.asarray(yv3))/luminosityNorm
    yv4=abs(np.asarray(yv4))/luminosityNorm
    yv5=abs(np.asarray(yv5))/luminosityNorm
    
    MAx1 = MovingAverage(xv1,10)
    MAx2 = MovingAverage(xv2,10)
    MAx5 = MovingAverage(xv5,10)
    MAy1 = MovingAverage(yv1,10)
    MAy2 = MovingAverage(yv2,10)
    MAy5 = MovingAverage(yv5,10)
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'$L/L_\mathrm{Edd}$',**font)
    plt.xlabel(r'$t\,[\times 10^4 GM/c^3]$',**font)
    #pylab.xlim([0,1])
    pylab.ylim([1e-5,1])
    plt.yscale('log')
    plot1, = plt.plot(MAx1, MAy1, "k--", dashes=[6,4], color = 'red', mec = 'red', markersize = 6, linewidth = 2)
    plot2, = plt.plot(MAx2, MAy2, "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 2)
    #plot3, = plt.plot(xv3, yv3, "k--", dashes=[6,4,2,4], color = 'blue', mec = 'blue', markersize = 6, linewidth = 2)
    #plot4, = plt.plot(xv4, yv4, "k--", dashes=[8,6,2,4,2,4], color = 'cyan', mec = 'cyan', markersize = 6, linewidth = 2)
    plot5, = plt.plot(MAx5, MAy5, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    plt.legend((r'S01E',r'S1E',r'S10E'), 'upper right', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../figures/luminosity_Tave.eps', format='eps', transparent='True')
    plt.show()
    
    freq1 = 10.**np.linspace(1.,2.5,10)
    
    A = np.fft.rfft(yv1)
    fft1 = 2*np.abs(A)**2*(xv1[len(xv1)-1]-xv1[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq1 = np.fft.rfftfreq(len(xv1), (xv1[1] - xv1[0])*1e4*timeUnit)
    A = np.fft.rfft(yv1a)
    fft1a = 2*np.abs(A)**2*(xv1a[len(xv1a)-1]-xv1a[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq1a = np.fft.rfftfreq(len(xv1a), (xv1a[1] - xv1a[0])*1e4*timeUnit)
    A = np.fft.rfft(yv2)
    fft2 = 2*np.abs(A)**2*(xv2[len(xv2)-1]-xv2[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq2 = np.fft.rfftfreq(len(xv2), (xv2[1] - xv2[0])*1e4*timeUnit)
    A = np.fft.rfft(yv3)
    fft3 = 2*np.abs(A)**2*(xv3[len(xv3)-1]-xv3[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq3 = np.fft.rfftfreq(len(xv3), (xv3[1] - xv3[0])*1e4*timeUnit)
    A = np.fft.rfft(yv4)
    fft4 = 2*np.abs(A)**2*(xv4[len(xv4)-1]-xv4[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq4 = np.fft.rfftfreq(len(xv4), (xv4[1] - xv4[0])*1e4*timeUnit)
    A = np.fft.rfft(yv5)
    fft5 = 2*np.abs(A)**2*(xv5[len(xv5)-1]-xv5[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq5 = np.fft.rfftfreq(len(xv5), (xv5[1] - xv5[0])*1e4*timeUnit)
    
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'Power [(rms/mean)$^2$/Hz]',**font)
    plt.xlabel(r'$\nu\,[Hz]$',**font)
    pylab.xlim([0.4,1000])
    #pylab.ylim([1e-6,1e0])
    plt.xscale('log')
    plt.yscale('log')
    plot1, = plt.plot(fftfreq1, fft1, "k--", color = 'red', mec = 'red', dashes=[6,4], markersize = 6, linewidth = 2)
    plot1a, = plt.plot(fftfreq1a, fft1a, "k--", color = 'red', mec = 'red', dashes=[6,4], markersize = 6, linewidth = 1)
    #plot6, = plt.plot(freq1,5.*freq1**(-1), 'k-', linewidth = 1)
    #plt.text(30., 0.2, r'$\propto \nu^{-1}$',**font)
    plt.legend((r'S01E',), 'lower left', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../shakura_01E_figs/PDS.eps', format='eps', transparent='True')
    plt.show()
    
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'Power [(rms/mean)$^2$/Hz]',**font)
    plt.xlabel(r'$\nu\,[Hz]$',**font)
    pylab.xlim([0.4,1000])
    pylab.ylim([1e-6,1e0])
    plt.xscale('log')
    plt.yscale('log')
    plot5, = plt.plot(fftfreq5, fft5, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    #plot6, = plt.plot(freq1,5.*freq1**(-1), 'k-', linewidth = 1)
    #plt.text(30., 0.2, r'$\propto \nu^{-1}$',**font)
    plt.legend((r'S10E',), 'lower left', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../shakura_10E_figs/PDS.eps', format='eps', transparent='True')
    plt.show()
    
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'Frequency*Power [(rms/mean)$^2$]',**font)
    plt.xlabel(r'$\nu\,[Hz]$',**font)
    pylab.xlim([0.4,1000])
    pylab.ylim([1e-10,1e2])
    plt.xscale('log')
    plt.yscale('log')
    plot1, = plt.plot(fftfreq1, fftfreq1*fft1, "k--", color = 'red', mec = 'red', dashes=[6,4], markersize = 6, linewidth = 2)
    plot1a, = plt.plot(fftfreq1a, fftfreq1a*fft1a, "k--", color = 'red', mec = 'red', dashes=[6,4], markersize = 6, linewidth = 1)
    #plot6, = plt.plot(freq1,5.*freq1**(-1), 'k-', linewidth = 1)
    #plt.text(30., 0.2, r'$\propto \nu^{-1}$',**font)
    plt.legend((r'S01E',), 'lower left', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../shakura_01E_figs/rms.eps', format='eps', transparent='True')
    plt.show()
    
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'Frequency*Power [(rms/mean)$^2$]',**font)
    plt.xlabel(r'$\nu\,[Hz]$',**font)
    pylab.xlim([0.4,1000])
    pylab.ylim([1e-4,1e0])
    plt.xscale('log')
    plt.yscale('log')
    plot5, = plt.plot(fftfreq5, fftfreq5*fft5, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    #plot6, = plt.plot(freq1,5.*freq1**(-1), 'k-', linewidth = 1)
    #plt.text(30., 0.2, r'$\propto \nu^{-1}$',**font)
    plt.legend((r'S10E',), 'lower left', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../shakura_10E_figs/rms.eps', format='eps', transparent='True')
    plt.show()
    return

#%% luminosity.py

def luminosity():
    """
    Inputs: none
    Return: none
    
    Takes in luminosity data and outputs two figures as 'luminosity.eps' and PDS.eps' respectivly
    """
    luminosityNorm = 2.29434e-21
    timeUnit = 3.26111e-05
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/luminosity','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/luminosity','r')]
    f3 = [np.array(line.split()).astype('float') for line in open('../shakura_3E_PP/luminosity','r')]
    f4 = [np.array(line.split()).astype('float') for line in open('../shakura_3Ep_PP/luminosity','r')]
    f5 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/luminosity','r')]
    
    xv1=[f1[i][0] for i in range(0,len(f1))]
    yv1=[f1[i][1] for i in range(0,len(f1))]
    
    xv2=[f2[i][0] for i in range(0,len(f2))]
    yv2=[f2[i][1] for i in range(0,len(f2))]
    
    xv3=[f3[i][0] for i in range(0,len(f3))]
    yv3=[f3[i][1] for i in range(0,len(f3))]
    
    xv4=[f4[i][0] for i in range(0,len(f4))]
    yv4=[f4[i][1] for i in range(0,len(f4))]
    
    xv5=[f5[i][0] for i in range(0,len(f5))]
    yv5=[f5[i][1] for i in range(0,len(f5))]
    
    xv1=np.asarray(xv1)/1e4
    xv2=np.asarray(xv2)/1e4
    xv3=np.asarray(xv3)/1e4
    xv4=np.asarray(xv4)/1e4
    xv5=np.asarray(xv5)/1e4
    
    yv1=abs(np.asarray(yv1))/luminosityNorm
    yv2=abs(np.asarray(yv2))/luminosityNorm
    yv3=abs(np.asarray(yv3))/luminosityNorm
    yv4=abs(np.asarray(yv4))/luminosityNorm
    yv5=abs(np.asarray(yv5))/luminosityNorm
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'$L/L_\mathrm{Edd}$',**font)
    plt.xlabel(r'$t\,[\times 10^4 GM/c^3]$',**font)
    #pylab.xlim([0,1])
    pylab.ylim([1e-5,1])
    plt.yscale('log')
    plot1, = plt.plot(xv1, yv1, "k--", dashes=[6,4], color = 'red', mec = 'red', markersize = 6, linewidth = 2)
    plot2, = plt.plot(xv2, yv2, "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 2)
    #plot3, = plt.plot(xv3, yv3, "k--", dashes=[6,4,2,4], color = 'blue', mec = 'blue', markersize = 6, linewidth = 2)
    #plot4, = plt.plot(xv4, yv4, "k--", dashes=[8,6,2,4,2,4], color = 'cyan', mec = 'cyan', markersize = 6, linewidth = 2)
    plot5, = plt.plot(xv5, yv5, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    plt.legend((r'S01E',r'S1E',r'S10E'), 'upper right', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../figures/luminosity.eps', format='eps', transparent='True')
    plt.show()
    
    freq1 = 10.**np.linspace(1.,2.5,10)
    
    A = np.fft.rfft(yv1)
    fft1 = 2*np.abs(A)**2*(xv1[len(xv1)-1]-xv1[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq1 = np.fft.rfftfreq(len(xv1), (xv1[1] - xv1[0])*1e4*timeUnit)
    A = np.fft.rfft(yv2)
    fft2 = 2*np.abs(A)**2*(xv2[len(xv2)-1]-xv2[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq2 = np.fft.rfftfreq(len(xv2), (xv2[1] - xv2[0])*1e4*timeUnit)
    A = np.fft.rfft(yv3)
    fft3 = 2*np.abs(A)**2*(xv3[len(xv3)-1]-xv3[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq3 = np.fft.rfftfreq(len(xv3), (xv3[1] - xv3[0])*1e4*timeUnit)
    A = np.fft.rfft(yv4)
    fft4 = 2*np.abs(A)**2*(xv4[len(xv4)-1]-xv4[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq4 = np.fft.rfftfreq(len(xv4), (xv4[1] - xv4[0])*1e4*timeUnit)
    A = np.fft.rfft(yv5)
    fft5 = 2*np.abs(A)**2*(xv5[len(xv5)-1]-xv5[0])*1e4*timeUnit/np.abs(A[0])**2
    fftfreq5 = np.fft.rfftfreq(len(xv5), (xv5[1] - xv5[0])*1e4*timeUnit)
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'Power [(rms/mean)$^2$/Hz]',**font)
    plt.xlabel(r'$\nu\,[Hz]$',**font)
    pylab.xlim([0.4,1000])
    pylab.ylim([1e-6,1e0])
    plt.xscale('log')
    plt.yscale('log')
    plot1, = plt.plot(fftfreq1, fft1, "k--", color = 'red', mec = 'red', dashes=[6,4], markersize = 6, linewidth = 2)
    plot2, = plt.plot(fftfreq2, fft2, "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 2)
    #plot3, = plt.plot(fftfreq3, fft3, "k--", dashes=[6,4,2,4], color = 'blue', mec = 'blue', markersize = 6, linewidth = 2)
    #plot4, = plt.plot(fftfreq4, fft4, "k--", dashes=[8,6,2,4,2,4], color = 'cyan', mec = 'cyan', markersize = 6, linewidth = 2)
    plot5, = plt.plot(fftfreq5, fft5, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    plot6, = plt.plot(freq1,5.*freq1**(-1), 'k-', linewidth = 1)
    plt.text(30., 0.2, r'$\propto \nu^{-1}$',**font)
    plt.legend((r'S01E',r'S1E',r'S10E'), 'lower left', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../figures/PDS.eps', format='eps', transparent='True')
    plt.show()
    return

#%% massFlux_calc.py
    
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

#%% massFluxXmin_LRBulkTest.py

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
#%% massFluxXmin_Tave.py

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

#%% massFluxXmin.py
    
def massFluxXmin():
    """
    Inputs: none
    Return: none
    
    Takes in mdot data and outputs figure labeled 'massFluxXmin.eps'
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
    plot1, = plt.plot(xv1, yv1, "k--", dashes=[6,4], color = 'red', mec = 'red', markersize = 6, linewidth = 2)
    plot2, = plt.plot(xv2, yv2, "k--", dashes=[2,3], color = 'green', mec = 'green', markersize = 6, linewidth = 2)
    plot3, = plt.plot(xv3, yv3, "k--", dashes=[10,8,2,4,2,4,2,4], color = 'orange', mec = 'orange', markersize = 6, linewidth = 2)
    plt.legend((r'S01E',r'S1E',r'S10E'), 'upper right', shadow=False, handlelength=3)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    plt.savefig('../figures/massFluxXmin.eps', format='eps', transparent='True')
    plt.show()
    return

#%% mdot.py

def mdot(file1,nt,nx,nskip,xmin,xmax,vmin,vmax,out):
    """
    Inputs: file1: filepath
            nt: int
            nx: int
            nskip: int 
            xmin: float
            xmax: float
            vmin: float
            vmax: float
            out: str (file path and name, e.g. ../figures/genericname)
    Return: none
    
    Takes in data files, returns and saves contour plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[9] correspond to inputs file1, ..., out.
    """
    if __name__ == "__main__":      #handels case where function called from command line
        mdot(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9])
        
    G      = 6.6725928e-8
    Msun   = 1.99e+33
    BHMass = 6.62
    c      = 2.99792458e+10
    
    munit  = BHMass*Msun
    lunit  = G*BHMass*Msun/c**2
    tunit  = lunit/c
    Medd   = 1.3*BHMass*1.e+17 #cgs
    
    f1 = open(file1,'r')
    alpha = 0.02
    
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
    
    zmin = 100
    zmax = 0
    for j in range(int(nt/nskip)):
        z.append([])
        for i in range(nx):
            if i == 0:
                y.append(yf[nskip*j*nx]/1.e4)
            if j == 0:
            	x.append(xf[i])
            z[j].append(zf[nskip*j*nx+i]*munit/tunit/Medd)
            if z[j][i] < zmin:
                zmin = z[j][i]
            if z[j][i] > zmax:
                zmax = z[j][i]
    
    # levels = MaxNLocator(nbins=100).bin_boundaries(zmin/5.,zmax/2.)
    levels = MaxNLocator(nbins=100).bin_boundaries(vmin,vmax)
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.xlabel(r'$R\,[GM/c^2]$',**font)                          #draws x-axis label
    plt.ylabel(r'$t\,[\times 10^4 GM/c^3]$',**font)             #draws y-axis label
    pylab.xlim([xmin,xmax])
    # pylab.ylim([0,y[nt-1]])
    pylab.ylim([2,3])
    # contours are *point* based plots, so convert our bound into point
    # centers
    plt.contourf(x, y, z, levels=levels, extend="both", cmap=cm.seismic)
    cbar = plt.colorbar(format="%.1f")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14)
    # cbar.set_label(r'$H/R$',**font)
    cbar.set_label(r'$\dot{m}$',**font)
    plt.savefig(out+".eps",format='eps',transparent='False')
    plt.show()
    return 

#%% Ptot0.py

def Ptot0(file1,nt,nx,nskip,xmin,xmax,out):
    """
    Inputs: file1: filepath
            nt: int
            nx: int
            nskip: int 
            xmin: float
            xmax: float
            out: str (file path and name, e.g. ../figures/genericname)
    Return: none
    
    Takes in data files, returns and saves contour plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[7] correspond to inputs file1, ..., out.
    """

    if __name__ == "__main__":       #handels case where function called from command line
        Ptot0(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])

    
    f1 = open(file1,'r')
    
    xf = []
    yf = []
    zf = []
    for line in f1:
        p1 = line.split()
        xf.append(float(p1[1]))
        yf.append(float(p1[0]))
        zf.append(np.log10(np.abs(float(p1[2]))))
       
    x = []
    y = []
    z = []
    
    zmin = 100
    zmax = -100
    for j in range(int(nt/nskip)):
        z.append([])
        for i in range(nx):
            if i == 0:
                y.append(yf[nskip*j*nx]/1e4)
            if j == 0:
                x.append(xf[i])
            z[j].append(zf[nskip*j*nx+i])
            if z[j][i] < zmin:
                zmin = z[j][i]
            if z[j][i] > zmax:
                zmax = z[j][i]
    
    levels = MaxNLocator(nbins=100).bin_boundaries(zmin,zmax)
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.xlabel(r'$R\,[GM/c^2]$',**font)                          #draws x-axis label
    plt.ylabel(r'$t\,[\times 10^4 GM/c^3]$',**font)             #draws y-axis label
    pylab.xlim([xmin,xmax])
    pylab.ylim([0,np.max(y)])
    # contours are *point* based plots, so convert our bound into point
    # centers
    plt.contourf(x, y, z, levels=levels, extend="both")
    cbar = plt.colorbar(format="%.2f")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14)
    cbar.set_label(r'$\log\,P_{z_0}$',**font)
    plt.savefig(out+".eps",format='eps',transparent='True')
    plt.show()
    return

#%% S-curve_open.py 

def Scurve_open():
    """
    Inputs: none
    Return: none
    
    Takes in Sigma and Tgas data, outputs and save figure
        """
    G=6.67e-8
    Msun=1.989e33
    c=2.99792485e10
    kb=1.3807e-16
    sigma=5.67e-5
    mp=1.67e-24
    aR=4.0*sigma/c
    
    alpha=0.02
    M=6.62*Msun
    bhspin = 0.
    mdot = [0.01,3.,10.]
    rg=G*M/(c**2)
    R=15.0*rg
    mu=0.615
    Kr=0.34
    
    S01nx=256
    if R==10.0*rg:
        S01istart=123
    else:
        S01istart=197 #r = 10, istart = 123; r = 15, istart = 197
    S01nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Tgas','r')]
    
    S01SigmaAvg = [0, 0, 0, 0, 0, 0]
    S01TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S01istart,S01istart+10):
        S01SigmaAvg[0] = S01SigmaAvg[0] + f1[i][2]/10
        S01TgasAvg[0]  = S01TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S01SigmaAvg)):
        for n in range(S01nstart,S01nstart+10):
            for i in range(S01istart,S01istart+10):
                S01SigmaAvg[a] = S01SigmaAvg[a] + f1[n*S01nx+i][2]/100
                S01TgasAvg[a]  = S01TgasAvg[a]  + f2[n*S01nx+i][2]/100
        S01nstart = S01nstart + 161
    
    S1nx=256
    if R==10.0*rg:
        S1istart=123
    else:
        S1istart=197 #r = 10, istart = 123; r = 15, istart = 197
    S1nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Tgas','r')]
    
    S1SigmaAvg = [0, 0, 0, 0, 0, 0]
    S1TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S1istart,S1istart+10):
        S1SigmaAvg[0] = S1SigmaAvg[0] + f1[i][2]/10
        S1TgasAvg[0]  = S1TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S1SigmaAvg)):
        for n in range(S1nstart,S1nstart+10):
            for i in range(S1istart,S1istart+10):
                S1SigmaAvg[a] = S1SigmaAvg[a] + f1[n*S1nx+i][2]/100
                S1TgasAvg[a]  = S1TgasAvg[a]  + f2[n*S1nx+i][2]/100
        S1nstart = S1nstart + 161
    
    S3nx=256
    if R==10.0*rg:
        S3istart=98
    else:
        S3istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S3nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_3E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_3E_PP/Tgas','r')]
    
    S3SigmaAvg = [0, 0, 0, 0, 0, 0]
    S3TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S3istart,S3istart+10):
        S3SigmaAvg[0] = S3SigmaAvg[0] + f1[i][2]/10
        S3TgasAvg[0]  = S3TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S3SigmaAvg)):
        for n in range(S3nstart,S3nstart+10):
            for i in range(S3istart,S3istart+10):
                S3SigmaAvg[a] = S3SigmaAvg[a] + f1[n*S3nx+i][2]/100
                S3TgasAvg[a]  = S3TgasAvg[a]  + f2[n*S3nx+i][2]/100
        S3nstart = S3nstart + 161
    
    S3nx = 256
    if R==10.0*rg:
        S3istart=98
    else:
        S3istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S3nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open(r'C:\\Users\\samps\\Desktop\\Py Scripts Project\\Data\\alpha', 'r')] #for line in open('../shakura_3Ep_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open(r'C:\\Users\\samps\\Desktop\\Py Scripts Project\\Data\\alpha', 'r')] #for line in open('../shakura_3Ep_PP/Tgas','r')]
    
    S3PSigmaAvg = [0, 0, 0, 0, 0, 0]
    S3PTgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S3istart,S3istart+10):
        S3PSigmaAvg[0] = S3PSigmaAvg[0] + f1[i][2]/10
        S3PTgasAvg[0]  = S3PTgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S3PSigmaAvg)):
        for n in range(S3nstart,S3nstart+10):
            for i in range(S3istart,S3istart+10):
                S3PSigmaAvg[a] = S3PSigmaAvg[a] + f1[n*S3nx+i][2]/100
                S3PTgasAvg[a]  = S3PTgasAvg[a]  + f2[n*S3nx+i][2]/100
        S3nstart = S3nstart + 161
    
    S10nx=256
    if R==10.0*rg:
        S10istart=98
    else:
        S10istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S10nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open(r'C:\\Users\\samps\\Desktop\\Py Scripts Project\\Data\\alpha', 'r')] #for line in open('../shakura_10E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open(r'C:\\Users\\samps\\Desktop\\Py Scripts Project\\Data\\alpha', 'r')] #for line in open('../shakura_10E_PP/Tgas','r')]
    
    S10SigmaAvg = [0, 0, 0, 0, 0, 0]
    S10TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S10istart,S10istart+10):
        S10SigmaAvg[0] = S10SigmaAvg[0] + f1[i][2]/10
        S10TgasAvg[0]  = S10TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S10SigmaAvg)):
        for n in range(S10nstart,S10nstart+10):
            for i in range(S10istart,S10istart+10):
                S10SigmaAvg[a] = S10SigmaAvg[a] + f1[n*S10nx+i][2]/100
                S10TgasAvg[a]  = S10TgasAvg[a]  + f2[n*S10nx+i][2]/100
        S10nstart = S10nstart + 161
    
    H = []
    SigmaDisk=[[3.25600e+03,3.25545e+03,3.25400e+03,3.27978e+03,4.33925e+03],
               [2.89744e+04,3.46920e+04,2.98701e+04,2.96050e+04],
               [8.63527e+03,4.85676e+03,8.48039e+03,1.16699e+04,1.37747e+04,6.66683e+03]]
    Tgas=[[1.62770e+06,2.86969e+06,2.92055e+06,4.16175e+06,4.26575e+06],
          [2.79669e+07,1.54046e+07,1.41224e+07,1.55943e+07],
          [2.86917e+07,1.06522e+07,9.38209e+06,9.22887e+06,1.12289e+07,1.11367e+07]]
    
    z1  = 1.0 + pow((1. - bhspin*bhspin), (1./3.))*(pow((1. + bhspin), (1./3.))+pow((1. - bhspin), (1./3.)))
    z2  = sqrt(3.*bhspin*bhspin + z1*z1)
    rms = (3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)));
    y = sqrt(R/rg)
    y0 = sqrt(rms)
    y1 = 2.*cos((acos(bhspin)-pi)/3.)
    y2 = 2.*cos((acos(bhspin)+pi)/3.)
    y3 = -2.*cos((acos(bhspin))/3.)
    m1 = M/Msun
    Ac = 1.+bhspin*bhspin*pow(y,-4)+2.*bhspin*bhspin*pow(y,-6)
    Bc = 1.+bhspin*pow(y,-3)
    Cc = 1.-3./y/y+2.*bhspin*pow(y,-3)
    Dc = 1.-2./y/y+bhspin*bhspin*pow(y,-4)
    Fc = 1.-2.*bhspin/y/y/y+bhspin*bhspin/y/y/y/y
    Gc = 1.-2./y/y+bhspin/y/y/y
    Rc = Fc*Fc/Cc-bhspin*bhspin/y/y*(Gc/sqrt(Cc)-1.)
    Ec = Ac*Ac/Bc/Bc*Cc/Dc*Rc
    Q0 = Bc/sqrt(Cc)/y
    Qc = Q0*(y-y0-1.5*bhspin*log(y/y0)-3.*pow(y1-bhspin,2)/y1/(y1-y2)/(y1-y3)*log((y-y1)/(y0-y1))-3.*pow(y2-bhspin,2)/y2/(y2-y1)/(y2-y3)*log((y-y2)/(y0-y2))-3.*pow(y3-bhspin,2)/y3/(y3-y1)/(y3-y2)*log((y-y3)/(y0-y3)))
    # for i in range(len(mdot)):
    #    if mdot[i] < 0.02:
    ## Middle region
    #        SigmaDisk.append(9.e4*pow(alpha,-0.8)*pow(m1,0.2)*pow(mdot[i],0.6)*pow(y,-1.2)*pow(Bc,-0.6)*sqrt(Cc)*pow(Dc,-0.8)*pow(Qc,0.6))
    #        Tgas.append(7.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(mdot[i],0.4)*pow(y,-1.8)*pow(Bc,-0.4)*pow(Dc,-0.2)*pow(Qc,0.4))
    #        H.append(1.e3*pow(alpha,-0.1)*pow(m1,0.9)*pow(mdot[i],0.2)*pow(y,2.1)*Ac*pow(Bc,-1.2)*sqrt(Cc)*pow(Dc,-0.6)*pow(Ec,-0.5)*pow(Qc,0.2))
    ## Outer region
    #        SigmaDisk.append(4.e5*pow(alpha,-0.8)*pow(m1,0.2)*pow(mdot[i],0.7)*pow(y,-1.5)*pow(Ac,0.1)*pow(Bc,-0.8)*sqrt(Cc)*pow(Dc,-0.85)*pow(Ec,-0.05)*pow(Qc,0.7))
    #        Tgas.append(2.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(mdot[i],0.3)*pow(y,-1.5)*pow(Ac,-0.1)*pow(Bc,-0.2)*pow(Dc,-0.15)*pow(Ec,0.05)*pow(Qc,0.3))
    #        H.append(4.e2*pow(alpha,-0.1)*pow(m1,0.9)*pow(mdot[i],0.15)*pow(y,2.25)*pow(Ac,0.95)*pow(Bc,-1.1)*sqrt(Cc)*pow(Dc,-0.575)*pow(Ec,-0.475)*pow(Qc,0.15))
    #    else:
    #        SigmaDisk.append(5./alpha/mdot[i]*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc)
    #        Tgas.append(5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-1.5)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    #        H.append(1.e5*m1*mdot[i]*Ac*Ac/Bc/Bc/Bc*sqrt(Cc)/Dc/Ec*Qc)
    #    print H[i]/rg
    print ("0.01E", 9.e4*pow(alpha,-0.8)*pow(m1,0.2)*pow(0.01,0.6)*pow(y,-1.2)*pow(Bc,-0.6)*sqrt(Cc)*pow(Dc,-0.8)*pow(Qc,0.6), 7.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(0.01,0.4)*pow(y,-1.8)*pow(Bc,-0.4)*pow(Dc,-0.2)*pow(Qc,0.4))
    print ("1E", 5./alpha/1.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    print ("3E", 5./alpha/3.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    print ("10E", 5./alpha/10.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    
    Tc=np.arange(3.0e5,4.19e7,1.0e2)
    omega=sqrt(G*M/(R**3))
    #T=Tc - (0.5*10**7)
    Sigma=np.sqrt(1.6*mp/kb*(32.*sigma/(27.*Kr*alpha*omega)*Tc**3-np.sqrt(32.*sigma/(27.*Kr*alpha*omega))/omega*4.*sigma/3./c*Tc**5))
    #S=Sigma +(1*10**4)
    
    ######################################
    # Novikov Thorne 1973 solutions
    ######################################
    Sigma_NT = np.arange(10.,1.e5,1.e2)
    Tin  = 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25)*Sigma_NT/Sigma_NT
    Tmid = 3.5e5*pow(alpha,1./3.)*pow(m1,-1./3.)/y*pow(Cc,-1./3.)*pow(Dc,1./3.)*pow(Sigma_NT,2./3.)
    ######################################
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel(r'$T_c\,\mathrm{[K]}$',**font)
    plt.xlabel(r'$\Sigma\,\mathrm{[g/cm^2]}$',**font)
    pylab.xlim([1200,1e5])
    pylab.ylim([1e6,1e8])
    plt.xscale('log')
    plt.yscale('log')
    plot2 = plt.plot(S01SigmaAvg[0], S01TgasAvg[0], "ko", color = 'red', mec = 'red', markersize = 6, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S1SigmaAvg[0], S1TgasAvg[0], "ks", color = 'green', mec = 'green', markersize = 6, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3SigmaAvg[0], S3TgasAvg[0], "kv", color = 'blue', mec = 'blue', markersize = 6, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3PSigmaAvg[0], S3PTgasAvg[0], "k<", color = 'cyan', mec = 'cyan', markersize = 6, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S10SigmaAvg[0], S10TgasAvg[0], "kD", color = 'orange', mec = 'orange', markersize = 6, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S01SigmaAvg[1], S01TgasAvg[1], "ko", color = 'red', mec = 'red', markersize = 8, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S01SigmaAvg[2], S01TgasAvg[2], "ko", color = 'red', mec = 'red', markersize = 10, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S01SigmaAvg[3], S01TgasAvg[3], "ko", color = 'red', mec = 'red', markersize = 12, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S01SigmaAvg[4], S01TgasAvg[4], "ko", color = 'red', mec = 'red', markersize = 14, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S01SigmaAvg[5], S01TgasAvg[5], "ko", color = 'red', mec = 'red', markersize = 16, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S1SigmaAvg[1], S1TgasAvg[1], "ks", color = 'green', mec = 'green', markersize = 8, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S1SigmaAvg[2], S1TgasAvg[2], "ks", color = 'green', mec = 'green', markersize = 10, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S1SigmaAvg[3], S1TgasAvg[3], "ks", color = 'green', mec = 'green', markersize = 12, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S1SigmaAvg[4], S1TgasAvg[4], "ks", color = 'green', mec = 'green', markersize = 14, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S1SigmaAvg[5], S1TgasAvg[5], "ks", color = 'green', mec = 'green', markersize = 16, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3SigmaAvg[1], S3TgasAvg[1], "kv", color = 'blue', mec = 'blue', markersize = 8, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3SigmaAvg[2], S3TgasAvg[2], "kv", color = 'blue', mec = 'blue', markersize = 10, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3SigmaAvg[3], S3TgasAvg[3], "kv", color = 'blue', mec = 'blue', markersize = 12, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3SigmaAvg[4], S3TgasAvg[4], "kv", color = 'blue', mec = 'blue', markersize = 14, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3SigmaAvg[5], S3TgasAvg[5], "kv", color = 'blue', mec = 'blue', markersize = 16, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3PSigmaAvg[1], S3PTgasAvg[1], "k<", color = 'cyan', mec = 'cyan', markersize = 8, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3PSigmaAvg[2], S3PTgasAvg[2], "k<", color = 'cyan', mec = 'cyan', markersize = 10, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3PSigmaAvg[3], S3PTgasAvg[3], "k<", color = 'cyan', mec = 'cyan', markersize = 12, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3PSigmaAvg[4], S3PTgasAvg[4], "k<", color = 'cyan', mec = 'cyan', markersize = 14, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S3PSigmaAvg[5], S3PTgasAvg[5], "k<", color = 'cyan', mec = 'cyan', markersize = 16, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S10SigmaAvg[1], S10TgasAvg[1], "kD", color = 'orange', mec = 'orange', markersize = 8, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S10SigmaAvg[2], S10TgasAvg[2], "kD", color = 'orange', mec = 'orange', markersize = 10, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S10SigmaAvg[3], S10TgasAvg[3], "kD", color = 'orange', mec = 'orange', markersize = 12, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S10SigmaAvg[4], S10TgasAvg[4], "kD", color = 'orange', mec = 'orange', markersize = 14, markerfacecolor = 'none', markeredgewidth = 1)
    plot2 = plt.plot(S10SigmaAvg[5], S10TgasAvg[5], "kD", color = 'orange', mec = 'orange', markersize = 16, markerfacecolor = 'none', markeredgewidth = 1)
    plot1 = plt.plot(Sigma, Tc, 'k', lw = 2)
    plot3 = plt.plot(Sigma_NT,Tin,'m--',lw = 2)
    plot4 = plt.plot(Sigma_NT,Tmid,'m--',lw = 2)
    
    plt.text(6.e3, 5.e7, r'$Q^+ > Q^-$',**font)
    plt.text(2.e3, 1.e7, r'$Q^+ < Q^-$',**font)
    plt.text(7.e3, 2.5e6, r'$Q^+ > Q^-$',**font)
    
    #plt.ylim(0,1e5)
    plt.legend((r'S01E',r'S1E',r'S3E',r'S3Ep',r'S10E'), 'lower right', shadow=False, numpoints = 1)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    if R==10.0*rg:
        plt.savefig('../figures/S-curve_10rg_open.eps', format='eps', transparent='True')
    else:
        plt.savefig('../figures/S-curve_15rg_open.eps', format='eps', transparent='True')
    plt.show()
    return 

#%% S-curve.py
def Scurve():
    """
    Inputs: none
    Return: none
    
    Takes in Sigma and Tgas data, outpus and save figure
    """
    G=6.67e-8
    Msun=1.989e33
    c=2.99792485e10
    kb=1.3807e-16
    sigma=5.67e-5
    mp=1.67e-24
    aR=4.0*sigma/c
    
    alpha=0.02
    M=6.62*Msun
    bhspin = 0.
    mdot = [0.01,3.,10.]
    rg=G*M/(c**2)
    R=10.0*rg
    mu=0.615
    Kr=0.34
    
    S01nx=256
    if R==10.0*rg:
        S01istart=123
    else:
        S01istart=197 #r = 10, istart = 123; r = 15, istart = 197
    S01nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_01E_PP/Tgas','r')]
    
    S01SigmaAvg = [0, 0, 0, 0, 0, 0]
    S01TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S01istart,S01istart+10):
        S01SigmaAvg[0] = S01SigmaAvg[0] + f1[i][2]/10
        S01TgasAvg[0]  = S01TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S01SigmaAvg)):
        for n in range(S01nstart,S01nstart+10):
            for i in range(S01istart,S01istart+10):
                S01SigmaAvg[a] = S01SigmaAvg[a] + f1[n*S01nx+i][2]/100
                S01TgasAvg[a]  = S01TgasAvg[a]  + f2[n*S01nx+i][2]/100
        S01nstart = S01nstart + 161
    
    S1nx=256
    if R==10.0*rg:
        S1istart=123
    else:
        S1istart=197 #r = 10, istart = 123; r = 15, istart = 197
    S1nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_1E_PP/Tgas','r')]
    
    S1SigmaAvg = [0, 0, 0, 0, 0, 0]
    S1TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S1istart,S1istart+10):
        S1SigmaAvg[0] = S1SigmaAvg[0] + f1[i][2]/10
        S1TgasAvg[0]  = S1TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S1SigmaAvg)):
        for n in range(S1nstart,S1nstart+10):
            for i in range(S1istart,S1istart+10):
                S1SigmaAvg[a] = S1SigmaAvg[a] + f1[n*S1nx+i][2]/100
                S1TgasAvg[a]  = S1TgasAvg[a]  + f2[n*S1nx+i][2]/100
        S1nstart = S1nstart + 161
    
    S3nx=256
    if R==10.0*rg:
        S3istart=98
    else:
        S3istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S3nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_3E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_3E_PP/Tgas','r')]
    
    S3SigmaAvg = [0, 0, 0, 0, 0, 0]
    S3TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S3istart,S3istart+10):
        S3SigmaAvg[0] = S3SigmaAvg[0] + f1[i][2]/10
        S3TgasAvg[0]  = S3TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S3SigmaAvg)):
        for n in range(S3nstart,S3nstart+10):
            for i in range(S3istart,S3istart+10):
                S3SigmaAvg[a] = S3SigmaAvg[a] + f1[n*S3nx+i][2]/100
                S3TgasAvg[a]  = S3TgasAvg[a]  + f2[n*S3nx+i][2]/100
        S3nstart = S3nstart + 161
    
    S3nx = 256
    if R==10.0*rg:
        S3istart=98
    else:
        S3istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S3nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_3Ep_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_3Ep_PP/Tgas','r')]
    
    S3PSigmaAvg = [0, 0, 0, 0, 0, 0]
    S3PTgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S3istart,S3istart+10):
        S3PSigmaAvg[0] = S3PSigmaAvg[0] + f1[i][2]/10
        S3PTgasAvg[0]  = S3PTgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S3PSigmaAvg)):
        for n in range(S3nstart,S3nstart+10):
            for i in range(S3istart,S3istart+10):
                S3PSigmaAvg[a] = S3PSigmaAvg[a] + f1[n*S3nx+i][2]/100
                S3PTgasAvg[a]  = S3PTgasAvg[a]  + f2[n*S3nx+i][2]/100
        S3nstart = S3nstart + 161
    
    S10nx=256
    if R==10.0*rg:
        S10istart=98
    else:
        S10istart=143 #r = 10, istart = 98; r = 15, istart = 143
    S10nstart=0
    
    f1 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/Sigma','r')]
    f2 = [np.array(line.split()).astype('float') for line in open('../shakura_10E_PP/Tgas','r')]
    
    S10SigmaAvg = [0, 0, 0, 0, 0, 0]
    S10TgasAvg  = [0, 0, 0, 0, 0, 0]
    for i in range(S10istart,S10istart+10):
        S10SigmaAvg[0] = S10SigmaAvg[0] + f1[i][2]/10
        S10TgasAvg[0]  = S10TgasAvg[0]  + f2[i][2]/10
    for a in range(1,len(S10SigmaAvg)):
        for n in range(S10nstart,S10nstart+10):
            for i in range(S10istart,S10istart+10):
                S10SigmaAvg[a] = S10SigmaAvg[a] + f1[n*S10nx+i][2]/100
                S10TgasAvg[a]  = S10TgasAvg[a]  + f2[n*S10nx+i][2]/100
        S10nstart = S10nstart + 161
    
    H = []
    SigmaDisk=[[3.25600e+03,3.25545e+03,3.25400e+03,3.27978e+03,4.33925e+03],
               [2.89744e+04,3.46920e+04,2.98701e+04,2.96050e+04],
               [8.63527e+03,4.85676e+03,8.48039e+03,1.16699e+04,1.37747e+04,6.66683e+03]]
    Tgas=[[1.62770e+06,2.86969e+06,2.92055e+06,4.16175e+06,4.26575e+06],
          [2.79669e+07,1.54046e+07,1.41224e+07,1.55943e+07],
          [2.86917e+07,1.06522e+07,9.38209e+06,9.22887e+06,1.12289e+07,1.11367e+07]]
    
    z1  = 1.0 + pow((1. - bhspin*bhspin), (1./3.))*(pow((1. + bhspin), (1./3.))+pow((1. - bhspin), (1./3.)))
    z2  = sqrt(3.*bhspin*bhspin + z1*z1)
    rms = (3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)));
    y = sqrt(R/rg)
    y0 = sqrt(rms)
    y1 = 2.*cos((acos(bhspin)-pi)/3.)
    y2 = 2.*cos((acos(bhspin)+pi)/3.)
    y3 = -2.*cos((acos(bhspin))/3.)
    m1 = M/Msun
    Ac = 1.+bhspin*bhspin*pow(y,-4)+2.*bhspin*bhspin*pow(y,-6)
    Bc = 1.+bhspin*pow(y,-3)
    Cc = 1.-3./y/y+2.*bhspin*pow(y,-3)
    Dc = 1.-2./y/y+bhspin*bhspin*pow(y,-4)
    Fc = 1.-2.*bhspin/y/y/y+bhspin*bhspin/y/y/y/y
    Gc = 1.-2./y/y+bhspin/y/y/y
    Rc = Fc*Fc/Cc-bhspin*bhspin/y/y*(Gc/sqrt(Cc)-1.)
    Ec = Ac*Ac/Bc/Bc*Cc/Dc*Rc
    Q0 = Bc/sqrt(Cc)/y
    Qc = Q0*(y-y0-1.5*bhspin*log(y/y0)-3.*pow(y1-bhspin,2)/y1/(y1-y2)/(y1-y3)*log((y-y1)/(y0-y1))-3.*pow(y2-bhspin,2)/y2/(y2-y1)/(y2-y3)*log((y-y2)/(y0-y2))-3.*pow(y3-bhspin,2)/y3/(y3-y1)/(y3-y2)*log((y-y3)/(y0-y3)))
    # for i in range(len(mdot)):
    #    if mdot[i] < 0.02:
    ## Middle region
    #        SigmaDisk.append(9.e4*pow(alpha,-0.8)*pow(m1,0.2)*pow(mdot[i],0.6)*pow(y,-1.2)*pow(Bc,-0.6)*sqrt(Cc)*pow(Dc,-0.8)*pow(Qc,0.6))
    #        Tgas.append(7.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(mdot[i],0.4)*pow(y,-1.8)*pow(Bc,-0.4)*pow(Dc,-0.2)*pow(Qc,0.4))
    #        H.append(1.e3*pow(alpha,-0.1)*pow(m1,0.9)*pow(mdot[i],0.2)*pow(y,2.1)*Ac*pow(Bc,-1.2)*sqrt(Cc)*pow(Dc,-0.6)*pow(Ec,-0.5)*pow(Qc,0.2))
    ## Outer region
    #        SigmaDisk.append(4.e5*pow(alpha,-0.8)*pow(m1,0.2)*pow(mdot[i],0.7)*pow(y,-1.5)*pow(Ac,0.1)*pow(Bc,-0.8)*sqrt(Cc)*pow(Dc,-0.85)*pow(Ec,-0.05)*pow(Qc,0.7))
    #        Tgas.append(2.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(mdot[i],0.3)*pow(y,-1.5)*pow(Ac,-0.1)*pow(Bc,-0.2)*pow(Dc,-0.15)*pow(Ec,0.05)*pow(Qc,0.3))
    #        H.append(4.e2*pow(alpha,-0.1)*pow(m1,0.9)*pow(mdot[i],0.15)*pow(y,2.25)*pow(Ac,0.95)*pow(Bc,-1.1)*sqrt(Cc)*pow(Dc,-0.575)*pow(Ec,-0.475)*pow(Qc,0.15))
    #    else:
    #        SigmaDisk.append(5./alpha/mdot[i]*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc)
    #        Tgas.append(5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-1.5)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    #        H.append(1.e5*m1*mdot[i]*Ac*Ac/Bc/Bc/Bc*sqrt(Cc)/Dc/Ec*Qc)
    #    print H[i]/rg
    print ("0.01E", 9.e4*pow(alpha,-0.8)*pow(m1,0.2)*pow(0.01,0.6)*pow(y,-1.2)*pow(Bc,-0.6)*sqrt(Cc)*pow(Dc,-0.8)*pow(Qc,0.6), 7.e8*pow(alpha,-0.2)*pow(m1,-0.2)*pow(0.01,0.4)*pow(y,-1.8)*pow(Bc,-0.4)*pow(Dc,-0.2)*pow(Qc,0.4))
    print ("1E", 5./alpha/1.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    print ("3E", 5./alpha/3.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    print ("10E", 5./alpha/10.*pow(y,3)/Ac/Ac*Bc*Bc*Bc*sqrt(Cc)*Ec/Qc, 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25))
    
    Tc=np.arange(3.0e5,4.19e7,1.0e2)
    omega=sqrt(G*M/(R**3))
    #T=Tc - (0.5*10**7)
    Sigma=np.sqrt(1.6*mp/kb*(32.*sigma/(27.*Kr*alpha*omega)*Tc**3-np.sqrt(32.*sigma/(27.*Kr*alpha*omega))/omega*4.*sigma/3./c*Tc**5))
    #S=Sigma +(1*10**4)
    
    ######################################
    # Novikov Thorne 1973 solutions
    ######################################
    Sigma_NT = np.arange(10.,1.e5,1.e2)
    Tin  = 5.e7*pow(alpha,-0.25)*pow(m1,-0.25)*pow(y,-0.75)/sqrt(Ac)*sqrt(Bc)*pow(Ec,0.25)*Sigma_NT/Sigma_NT
    Tmid = 3.5e5*pow(alpha,1./3.)*pow(m1,-1./3.)/y*pow(Cc,-1./3.)*pow(Dc,1./3.)*pow(Sigma_NT,2./3.)
    ######################################
    
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.tick_params(length=8)
    plt.tick_params(length=5,which='minor')
    plt.ylabel('$T_c\,\mathrm{[K]}$',**font)
    plt.xlabel('$\Sigma\,\mathrm{[g/cm^2]}$',**font)
    pylab.xlim([1200,1e5])
    pylab.ylim([1e6,1e8])
    plt.xscale('log')
    plt.yscale('log')
    plot2 = plt.plot(S01SigmaAvg[0], S01TgasAvg[0], "ko", color = 'red', mec = 'red', markersize = 6)
    plot2 = plt.plot(S1SigmaAvg[0], S1TgasAvg[0], "ks", color = 'green', mec = 'green', markersize = 6)
    plot2 = plt.plot(S3SigmaAvg[0], S3TgasAvg[0], "kv", color = 'blue', mec = 'blue', markersize = 6)
    plot2 = plt.plot(S3PSigmaAvg[0], S3PTgasAvg[0], "k<", color = 'cyan', mec = 'cyan', markersize = 6)
    plot2 = plt.plot(S10SigmaAvg[0], S10TgasAvg[0], "kD", color = 'orange', mec = 'orange', markersize = 6)
    plot2 = plt.plot(S01SigmaAvg[1], S01TgasAvg[1], "ko", color = 'red', mec = 'red', markersize = 8)
    plot2 = plt.plot(S01SigmaAvg[2], S01TgasAvg[2], "ko", color = 'red', mec = 'red', markersize = 10)
    plot2 = plt.plot(S01SigmaAvg[3], S01TgasAvg[3], "ko", color = 'red', mec = 'red', markersize = 12)
    plot2 = plt.plot(S01SigmaAvg[4], S01TgasAvg[4], "ko", color = 'red', mec = 'red', markersize = 14)
    plot2 = plt.plot(S01SigmaAvg[5], S01TgasAvg[5], "ko", color = 'red', mec = 'red', markersize = 16)
    plot2 = plt.plot(S1SigmaAvg[1], S1TgasAvg[1], "ks", color = 'green', mec = 'green', markersize = 8)
    plot2 = plt.plot(S1SigmaAvg[2], S1TgasAvg[2], "ks", color = 'green', mec = 'green', markersize = 10)
    plot2 = plt.plot(S1SigmaAvg[3], S1TgasAvg[3], "ks", color = 'green', mec = 'green', markersize = 12)
    plot2 = plt.plot(S1SigmaAvg[4], S1TgasAvg[4], "ks", color = 'green', mec = 'green', markersize = 14)
    plot2 = plt.plot(S1SigmaAvg[5], S1TgasAvg[5], "ks", color = 'green', mec = 'green', markersize = 16)
    plot2 = plt.plot(S3SigmaAvg[1], S3TgasAvg[1], "kv", color = 'blue', mec = 'blue', markersize = 8)
    plot2 = plt.plot(S3SigmaAvg[2], S3TgasAvg[2], "kv", color = 'blue', mec = 'blue', markersize = 10)
    plot2 = plt.plot(S3SigmaAvg[3], S3TgasAvg[3], "kv", color = 'blue', mec = 'blue', markersize = 12)
    plot2 = plt.plot(S3SigmaAvg[4], S3TgasAvg[4], "kv", color = 'blue', mec = 'blue', markersize = 14)
    plot2 = plt.plot(S3SigmaAvg[5], S3TgasAvg[5], "kv", color = 'blue', mec = 'blue', markersize = 16)
    plot2 = plt.plot(S3PSigmaAvg[1], S3PTgasAvg[1], "k<", color = 'cyan', mec = 'cyan', markersize = 8)
    plot2 = plt.plot(S3PSigmaAvg[2], S3PTgasAvg[2], "k<", color = 'cyan', mec = 'cyan', markersize = 10)
    plot2 = plt.plot(S3PSigmaAvg[3], S3PTgasAvg[3], "k<", color = 'cyan', mec = 'cyan', markersize = 12)
    plot2 = plt.plot(S3PSigmaAvg[4], S3PTgasAvg[4], "k<", color = 'cyan', mec = 'cyan', markersize = 14)
    plot2 = plt.plot(S3PSigmaAvg[5], S3PTgasAvg[5], "k<", color = 'cyan', mec = 'cyan', markersize = 16)
    plot2 = plt.plot(S10SigmaAvg[1], S10TgasAvg[1], "kD", color = 'orange', mec = 'orange', markersize = 8)
    plot2 = plt.plot(S10SigmaAvg[2], S10TgasAvg[2], "kD", color = 'orange', mec = 'orange', markersize = 10)
    plot2 = plt.plot(S10SigmaAvg[3], S10TgasAvg[3], "kD", color = 'orange', mec = 'orange', markersize = 12)
    plot2 = plt.plot(S10SigmaAvg[4], S10TgasAvg[4], "kD", color = 'orange', mec = 'orange', markersize = 14)
    plot2 = plt.plot(S10SigmaAvg[5], S10TgasAvg[5], "kD", color = 'orange', mec = 'orange', markersize = 16)
    plot1 = plt.plot(Sigma, Tc, 'k', lw = 2)
    plot3 = plt.plot(Sigma_NT,Tin,'m--',lw = 2)
    plot4 = plt.plot(Sigma_NT,Tmid,'m--',lw = 2)
    
    plt.text(6.e3, 5.e7, r'$Q^+ > Q^-$',**font)
    plt.text(2.e3, 1.e7, r'$Q^+ < Q^-$',**font)
    plt.text(7.e3, 2.5e6, r'$Q^+ > Q^-$',**font)
    
    #plt.ylim(0,1e5)
    plt.legend((r'S01E',r'S1E',r'S3E',r'S3Ep',r'S10E'), 'lower right', shadow=False, numpoints = 1)
    leg =  plt.gca().get_legend()
    leg.get_frame().set_alpha(0)
    ltext = leg.get_texts()
    plt.setp(ltext,**font)
    
    if R==10.0*rg:
        plt.savefig('../figures/S-curve_10rg.eps', format='eps', transparent='True')
    else:
        plt.savefig('../figures/S-curve_15rg.eps', format='eps', transparent='True')
    plt.show()
    return

#%% sigma_prad_t.py

def sigma_prad_t(file1, file2, nx, tstart, ymax1,ymin2, ymax2, out):
    """
    Inputs:
        file1: filepath
        file2: filepath
        nx: int
        tstart: int
        ymax1: float
        ymin2: float
        ymax2: float
        out: filepath
        
    Return: none
    
    Takes in data files, returns and saves plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[8] correspond to inputs file1, ..., out.
    """

    if __name__ == "__main__":       #handels case where function called from command line
        sigma_prad_t(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7], sys.argv[8])

    
    istart=tstart*nx
    istop=(tstart+1)*nx
    
    f1 = [np.array(line.split()).astype('float') for line in open(file1,'r')]
    f2 = [np.array(line.split()).astype('float') for line in open(file2,'r')]
    
    xv1=[f1[i][1] for i in range(istart,istop)]
    yv1=[f1[i][2] for i in range(istart,istop)]
    
    xv2=[f2[i][1] for i in range(istart,istop)]
    yv2=[np.log10(f2[i][2]) for i in range(istart,istop)]
    
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
    ax1.set_ylabel(r'$\log\,(P_\mathrm{rad}/P_\mathrm{gas})$',**font)
    ax1.set_ylim([ymin2,ymax2])
    ax1.set_xlim([xv2[0],25])
    #ax1.set_yscale('log')
    ax1.tick_params(length=5, which='minor')
    ax1.tick_params(length=8, which='major')
    axtext = ax1.yaxis.get_offset_text()
    plt.setp(axtext,**font)
    ax1.plot(xv2, yv2, "k--", dashes=[2,3], color = 'blue', mec = 'blue', markersize = 6, linewidth = 2)
    
    
    #ax0.legend((r'288x96x96_2level',r'384x128x128_2level',r'288x96x96_3level',r'384x128x128_3level'), 'upper left', shadow=False, handlelength=3)
    #leg = ax0.get_legend()
    #ltext  = leg.get_texts()
    #plt.setp(ltext, **font2)
    axtext = ax0.get_xmajorticklabels()
    plt.setp(axtext,**font)
    axtext = ax0.get_ymajorticklabels()
    plt.setp(axtext,**font)
    axtext = ax1.get_ymajorticklabels()
    plt.setp(axtext,**font)
    
    plt.savefig(out+".eps", format='eps', transparent='True')
    plt.show()
    return

#%% sigma_stress_t.py
    
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

#%% sigma_tau_t.py

def sigma_tau_t(file1, file2, nx, tstart, ymax1, ymax2, out):
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

    if __name__ == "__main__":       #handels case where function called from command line
        sigma_tau_t(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7])

    
    istart=tstart*nx
    istop=(tstart+1)*nx
    
    f1 = [np.array(line.split()).astype('float') for line in open(file1,'r')]
    f2 = [np.array(line.split()).astype('float') for line in open(file2,'r')]
    
    xv1=[f1[i][1] for i in range(istart,istop)]
    yv1=[f1[i][2] for i in range(istart,istop)]
    
    xv2=[f2[i][1] for i in range(istart,istop)]
    yv2=[f2[i][2] for i in range(istart,istop)]
    
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
    ax1.set_ylabel(r'$\tau$',**font)
    ax1.set_ylim([0,ymax2])
    ax1.set_xlim([xv2[0],25])
    #ax1.set_yscale('log')
    ax1.tick_params(length=5, which='minor')
    ax1.tick_params(length=8, which='major')
    axtext = ax1.yaxis.get_offset_text()
    plt.setp(axtext,**font)
    ax1.plot(xv2, yv2, "k--", dashes=[2,3], color = 'blue', mec = 'blue', markersize = 6, linewidth = 2)
    
    
    #ax0.legend((r'288x96x96_2level',r'384x128x128_2level',r'288x96x96_3level',r'384x128x128_3level'), 'upper left', shadow=False, handlelength=3)
    #leg = ax0.get_legend()
    #ltext  = leg.get_texts()
    #plt.setp(ltext, **font2)
    axtext = ax0.get_xmajorticklabels()
    plt.setp(axtext,**font)
    axtext = ax0.get_ymajorticklabels()
    plt.setp(axtext,**font)
    axtext = ax1.get_ymajorticklabels()
    plt.setp(axtext,**font)
    
    plt.savefig(out+".eps", format='eps', transparent='True')
    plt.show()
    return

#%% tau.py

def tau(file1, nt, nx, nskip, xmin,xmax, vmax, out):
    """
    Inputs:
        file1: filepath
        nt: int
        nx: int
        nskip: int
        xmin: float
        xmax: float
        vmax: float
        out: filepath
        
    Return: none
    
    Takes in data files, returns and saves plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[8] correspond to inputs file1, ..., out.
    """

    if __name__ == "__main__":       #handels case where function called from command line
        tau(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7], sys.argv[8])


    f1 = open(file1,'r')
    
    alpha = 0.02
    pi    = np.pi
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
            z[j].append(zf[nskip*j*nx+i])
            if z[j][i] < zmin:
                zmin = z[j][i]
            if z[j][i] > zmax:
                zmax = z[j][i]
    
    levels = MaxNLocator(nbins=100).bin_boundaries(zmin,vmax)
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
    cbar = plt.colorbar(format="%.0f")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14)
    cbar.set_label(r'$\tau$',**font)
    plt.savefig(out+".eps",format='eps',transparent='True')
    plt.show()
    return
 
#%% tgas.py

def tgas(file1, nt, nx, nskip, xmin,xmax, vmax, out):
    """
    Inputs:
        file1: filepath
        nt: int
        nx: int
        nskip: int
        xmin: float
        xmax: float
        vmax: float
        out: filepath
        
    Return: none
    
    Takes in data files, returns and saves plot as out+'.eps'
    If called from command line, sys.argv[1], ..., sys.argv[8] correspond to inputs file1, ..., out.
    """

    if __name__ == "__main__":       #handels case where function called from command line
        tgas(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7], sys.argv[8])

    
    f1 = open(file1,'r')
    xf = []
    yf = []
    zf = []
    for line in f1:
        p1 = line.split()
        xf.append(float(p1[1]))
        yf.append(float(p1[0]))
        zf.append(np.abs(float(p1[2])))
       
    x = []
    y = []
    z = []
    
    zmin = 100
    zmax = -100
    for j in range(int(nt/nskip)):
        z.append([])
        for i in range(nx):
            if i == 0:
                y.append(yf[nskip*j*nx]/1.e4)
            if j == 0:
                x.append(xf[i])
            z[j].append(zf[nskip*j*nx+i])
            if z[j][i] < zmin:
                zmin = z[j][i]
            if z[j][i] > zmax:
                zmax = z[j][i]
    
    levels = MaxNLocator(nbins=100).bin_boundaries(zmin,vmax)
    font = {'fontname':'Times New Roman','fontsize':16, 'weight':'normal'}
    plt.xticks(**font)
    plt.yticks(**font)
    plt.xlabel(r'$R\,[GM/c^2]$',**font)                          #draws x-axis label
    plt.ylabel(r'$t\,[GM/c^3]$',**font)             #draws y-axis label
    pylab.xlim([xmin,xmax])
    # pylab.ylim([0,np.max(y)])
    # contours are *point* based plots, so convert our bound into point
    # centers
    plt.contourf(x, y, z, levels=levels, extend="both")
    cbar = plt.colorbar(format="%.1e")
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family("Times New Roman")
        l.set_size(14)
    cbar.set_label(r'$T_\mathrm{gas}\,[\mathrm{K}]$',**font)
    # cbar.set_label(r'$\log(T_\mathrm{rad})\,[\mathrm{K}]$',**font)
    plt.savefig(out+".eps",format='eps',transparent='True')
    plt.show()
    return

#%% vertStress.py
    
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

#%% zdensity.py

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