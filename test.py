from __future__ import division
from scipy.integrate import quad
from math import *
import pyfits
import glob
import pylab as P
import numpy as N 
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches
from scipy.optimize import fmin
from operator import add
from numpy.linalg import inv
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
import sys


def plot_mwd(RA,Dec,org=0,title='Mollweide projection', projection='mollweide'):
    ''' RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    title is the title of the figure.
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    x = N.remainder(RA+360-org,360) # shift RA values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    tick_labels = N.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = N.remainder(tick_labels+360+org,360)
    fig = P.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=projection,)# axisbg ='LightCyan')
    ax.scatter(N.radians(x),N.radians(Dec))  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)






coord = SkyCoord([2,0],[0,2], frame='icrs', unit='deg', )
	
gc=coord.galactic

gcl=gc.l.degree
gcb=gc.b.degree

ra = coord.ra.degree
dec = coord.dec.degree

print type(ra)
print dec
sys.exit()

P.figure()
P.scatter(ra,dec)
P.grid()

coord = N.array([(0,30), (60,-45), (240,15), (150,-75)])
plot_mwd(coord[:,0],coord[:,1], org=90, title ='Test plot_mwd')

plot_mwd(ra,dec,org=0, title = 'oui?')



P.figure()
ax = P.subplot(111,projection="mollweide")
ax.scatter(ra,dec)
P.grid()


P.show()
