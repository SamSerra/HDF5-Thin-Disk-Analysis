import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab
"""
Input: none
Return: none

Takes in Mdot and luminosity data, manipulates and plots them, and saves the figure as efficiency.eps
"""
luminosityNorm = 2.29434e-21
massFluxNorm = 2.29434e-21
timeUnit = 3.26111e-05

f1 = [np.array(line.split()).astype('float') for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] #for line in open('../shakura_01E_PP/luminosity','r')]
f2 = [np.array(line.split()).astype('float') for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] #for line in open('../shakura_01E_PP/Mdot','r')]
f3 = [np.array(line.split()).astype('float') for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] #for line in open('../shakura_1E_PP/luminosity','r')]
f4 = [np.array(line.split()).astype('float') for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] #for line in open('../shakura_1E_PP/Mdot','r')]
f5 = [np.array(line.split()).astype('float') for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] #for line in open('../shakura_10E_PP/luminosity','r')]
f6 = [np.array(line.split()).astype('float') for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] #for line in open('../shakura_10E_PP/Mdot','r')]

xv1=[f1[i][0] for i in range(0,len(f1))]
yv1=[f1[i][1] for i in range(0,len(f1))]

xv2=[f2[i][0] for i in range(0,len(f2))] #for i in range(2,len(f2),256)]
yv2=[f2[i][2] for i in range(0,len(f2))] #for i in range(2,len(f2),256)]

xv3=[f3[i][0] for i in range(0,len(f3))]
yv3=[f3[i][1] for i in range(0,len(f3))]

xv4=[f4[i][0] for i in range(0,len(f4))] #for i in range(2,len(f4),256)]
yv4=[f4[i][2] for i in range(0,len(f4))] #for i in range(2,len(f4),256)]

xv5=[f5[i][0] for i in range(0,len(f5))]
yv5=[f5[i][1] for i in range(0,len(f5))]

xv6=[f6[i][0] for i in range(0,len(f6))] #for i in range(2,len(f6),256)]
yv6=[f6[i][2] for i in range(0,len(f6))] #for i in range(2,len(f6),256)]

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
plt.ylabel(u'$\eta$',**font)
plt.xlabel(u'$t\,[\times 10^4 GM/c^3]$',**font)
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
