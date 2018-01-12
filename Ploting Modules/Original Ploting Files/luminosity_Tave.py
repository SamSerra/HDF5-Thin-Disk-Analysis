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
