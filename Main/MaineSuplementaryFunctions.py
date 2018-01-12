import numpy as np
import matplotlib.pyplot as plt
"""
Module to Supplement Maine.py
"""

def ConvertToPolar(xArray,yArray):
	"""
	Takes in two arrays coresponding to x and y coordinates and transforms
	them to polar coordinates.

	Returns radius array, angle array
	"""
	assert xArray.shape == yArray.shape, "Array's must have same shape!"

	rad, theta = np.ones(xArray.shape), np.ones(xArray.shape)

	rad=np.sqrt(xArray**2+yArray**2)
	theta=np.arctan2(yArray,xArray)

	return rad, theta


def plot_polar_contour(values, theta, radius):
	"""
	Takes in z,theta,r arrays
	Generates a polar contour plot symetric about theta = 0
	"""
	theta = np.array(theta)
	radius = np.array(radius)

	values = np.array(values)
	values = values.reshape(len(theta), len(radius))

	r, tet = np.meshgrid(radius, theta)	#meshgrid set up opposite of values b/c 
										#pyplot uses transpose 

	print("Generating Plot...")	
	maxbound = np.max(theta)	#bounds for plot

	fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
	cax = ax.contourf(tet, r, values, 30, cmap = "hot", vmin = 1.0e-29, vmax = 1.0e-18)

	cb = fig.colorbar(cax)
	cb.set_label("Mass Density")

	x=ax.axes.get_xaxis()
	x.set_visible(False)

	plt.savefig('MassDensityContour.png', dpi = 600)
	print("Plot Saved")
	return fig, ax, cax