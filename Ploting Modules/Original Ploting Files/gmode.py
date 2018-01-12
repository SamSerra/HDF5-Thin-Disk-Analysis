import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
import pylab
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
"""
Input: none
Return: none

Takes in Mdot data, plots PDSMdot figure, and saves it as PDSmdot01E_r7p7_8p0_8p3.eps
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

f1 = [np.array(line.split()).astype('float') for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] 
#for line in open('../shakura_01E_PP/Mdot','r')]
f2 = [np.array(line.split()).astype('float')  for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] 
#for line in open('../shakura_1E_PP/Mdot','r')]
f3 = [np.array(line.split()).astype('float')  for line in open(r'C:\Users\samps\Desktop\Py Scripts Project\Data\alpha','r')] 
#for line in open('../shakura_10E_PP/Mdot','r')]


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
plt.xlabel(u'$\\nu\,[\\times 10^{-3}c^3/GM]$',**font)                          #draws x-axis label
plt.ylabel(u'$\mathrm{PDS}\,\dot{m}$',**font)
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