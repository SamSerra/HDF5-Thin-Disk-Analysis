import numpy as np
import matplotlib.pyplot as plt
import h5py

"""
Plots original grid and centers of each cell
"""

nDims = 2	#number of dimesions
nNodes = 1	#number of nodes per cell (nodes are verticies of cells)
if nDims == 1: nNodes = 2
elif nDims == 2: nNodes = 4
else: nNodes = 8

nSkip = 1


## File Addresses ##
fileExtension = r'\Users\samps\Desktop\HDF5 Project\Test Data\output\proc0000'
datasetName = 'Cosmos++'

## Data Handling ##
##---------------##

dataf = h5py.File(fileExtension + '\outHDF5-cycle0000', 'r')
data = dataf[datasetName]

## Get Chunk Data ##
chunk_shape = data.chunks 
chunk_size = chunk_shape[0]	


## Grid Handling ##
##---------------##

gridf = h5py.File(fileExtension + '\gridoutHDF5-cycle0000', 'r')
grid = gridf[datasetName]


## Read x,y Coords ##
xNodeCoords, yNodeCoords = np.ones(nNodes*chunk_size), np.ones(nNodes*chunk_size)

for n in np.arange(0,chunk_size*nNodes,nSkip):
	xNodeCoords[n] = grid[2*n]
	yNodeCoords[n] = grid[2*n+1]
	print('Reading Data from Node {} \n\tx-coord: {} \n\ty-coord: {}'\
		.format(n,xNodeCoords[n], yNodeCoords[n]))

## Compute Center of Cells ##
xAvg, yAvg = np.ones(chunk_size), np.ones(chunk_size)

ncnt = 0
for n in np.arange(chunk_size):
	xAvg[n] = 0.25*(xNodeCoords[ncnt]+xNodeCoords[ncnt+1]\
					+xNodeCoords[ncnt+2]+xNodeCoords[ncnt+3])
	yAvg[n] = 0.25*(yNodeCoords[ncnt]+yNodeCoords[ncnt+1]\
					+yNodeCoords[ncnt+2]+yNodeCoords[ncnt+3])
	ncnt += nNodes
	print('Center for cell {}: \n\t ({},{})'.format(n,xAvg[n],yAvg[n]))


f = plt.figure()
ax = f.gca()

ax.set_xlabel('Radius'+ r'[\frac {GM}{c^2}]')
ax.set_ylabel('Height'+ r'[\frac {GM}{c^2}]')
ax.set_title('Initial Grid (ploting every {} points)'.format(nSkip))

plt.scatter(xNodeCoords,yNodeCoords, marker = '.')
plt.scatter(xAvg,yAvg, c = 'g', marker = '.')

plt.show()
#plt.savefig('Initial Grid Plot.png')

dataf.close()