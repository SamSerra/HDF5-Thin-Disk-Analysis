import sys
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pylab
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
"""
Plots 
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
plt.xlabel(u'$r\,[GM/c^2]$',**font)                          #draws x-axis label
plt.ylabel(u'$\\kappa\,[\\times 10^{-3}\,c^3/GM]$',**font)             #draws y-axis label
#plt.xscale('log')
pylab.xlim([7.6,8.4])
pylab.ylim([3.47,3.52]) #RADP
#plt.savefig('../figures/gmode_cavity.eps',format='eps',Transparent='True')
plt.show()