import numpy as np
import matplotlib.pyplot as plt

"""
Calculates the distance from the origin to the center of each
cell
"""

f=np.loadtxt("massDensityTest.txt")
print("Data has {} Rows".format(np.size(f,0)))
rad = np.ones((np.size(f,0),))

for n in np.arange(np.size(f,0)):
	rad[n]=np.sqrt(f[n,0]**2+f[n,1]**2)
	print (rad[n])

with open('ShellMassDensity Test.txt', 'w') as f1:
	for i in np.arange(np.size(f,0)):
		f1.write('{}\t{}\t{}\n'.format(str(rad[i]),str(f[i,2]),str(rad[i]-rad[i-1])))
