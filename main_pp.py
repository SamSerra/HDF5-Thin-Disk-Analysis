"""
Reads Mass density and x,y coodinate from HDF5 data file
Note: reads as single-pression float
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import MaineSuplementaryFunctions

## Setup Matching Original Simulation ##
##------------------------------------##

nx=256	#256 divisions of radius, 192 divisions of theta, 1 division of phi
ny=192
nz=1

nDims = 2	#number of dimesions
nNodes = 1	#number of nodes per cell (nodes are verticies of cells)

if nDims == 1: nNodes = 2
elif nDims == 2: nNodes = 4
else: nNodes = 8


## Data Handling ##
##---------------##

## Open Data File and Dataset ##
fileExtension = r'\Users\samps\Projects\HDF5 Project\Test Data\output\proc0000'
datasetName = 'Cosmos++'

dataf = h5py.File(fileExtension + '\outHDF5-cycle0000', 'r')
data = dataf[datasetName]

## Get Chunk Data ##
chunk_shape = data.chunks 
chunk_size = chunk_shape[0]	

## Get Mass Density Data and Store in Array ##
mass_density = np.ones(chunk_shape) 
mass_density = data[chunk_size*2:chunk_size*3] #Mass density is the 3rd chunk


## Grid Handling ##
##---------------##

gridf = h5py.File(fileExtension + '\gridoutHDF5-cycle0000', 'r')
grid = gridf[datasetName]

## Read x,y Coords ##
xNodeCoords, yNodeCoords = np.ones(nNodes*chunk_size), np.ones(nNodes*chunk_size)

for n in np.arange(chunk_size*nNodes):
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

## Convert from Cartisian to Polar Coordinates ##
rad, theta = np.ones(len(xNodeCoords)), np.ones(len(xNodeCoords))
rad, theta = MaineSuplementaryFunctions.ConvertToPolar(xAvg,yAvg)


for line in np.arange(len(rad)):
		print("Radius: {}\t Theta: {}".format(rad[line],theta[line]))


## Write to Test File ##
##--------------------##
"""
with open('massDensityTest.txt', 'w') as f1:
	for line in np.arange(len(mass_density)):
		f1.write(str(xAvg[line]) +'\t' + str(yAvg[line]) + '\t' + \
			str(rad[line]) +'\t' + str(theta[line]) + '\t' + \
			str(mass_density[line]) + '\n')

"""
## Plot Mass Density on Polar Plot ##
##---------------------------------##

## Calculate Average Radius Per Shell ##

radAvg = np.ones(nx)

for i in range(nx):
	radAvg[i]=sum(rad[i::nx])/ny
	print("Average Radius in Shell {} is {}".format(i,radAvg[i]))

## Calculate Average Theta per Slice ##

thetaAvg = np.ones(ny)

ncnt = 0
for i in range(ny):
	thetaAvg[i] = sum(theta[ncnt:ncnt+nx])/nx
	ncnt += nx
	print("Average Theta for Line {} is {}".format(i,thetaAvg[i]))

## Put Mass Density in ny by nx Array ##

mass_density_array = np.ones((nx,ny))
ncnt = 0
for i in range(ny):
	for j in range(nx):
		mass_density_array[j][i] = mass_density[ncnt+j]
	ncnt += nx

mass_density_array = np.transpose(mass_density_array)

## Generate Polar Plot ##

MaineSuplementaryFunctions.plot_polar_contour(mass_density_array,thetaAvg,radAvg)

## Clean Up ## 
##----------##
dataf.close()
gridf.close()

