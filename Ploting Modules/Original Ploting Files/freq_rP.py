import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
import pylab
from scipy import signal
from matplotlib import rcParams
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import colors, ticker, cm
"""
Inputs: none
Return: none

Imports Mdot data, preforms DFT, and plots figure
"""
rcParams.update({'figure.autolayout': True})

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
plt.xlabel(u'$R\,[GM/c^2]$',**font)                          #draws x-axis label
plt.ylabel(u'$\nu\,[\times 10^{-3}c^3/GM]$',**font)
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
