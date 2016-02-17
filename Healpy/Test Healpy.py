import matplotlib.pyplot as P
import numpy as N 
import healpy as hp
import sys


nside = 2**5 #nsize must be 2**k. numbers of pixels = 12*(nsize**2)


m = N.arange(hp.nside2npix(nside))
hp.mollview(m,nest=None,title='first')

#Procedure to read and write maps form fits files
'''
hp.write_map("my_map.fits", m)
wmap_map_I = hp.read_map('my_map.fits')
hp.mollview(wmap_map_I,title='test')
'''


hp.mollview(m,title='Histogram equalized Ecliptic', unit='mK', xsize=2000)
hp.graticule() # print the meridians and parralels


#ca je ne comprend pas trop mais je ne risaue pas de m'en servir
"""
hp.gnomview(m, rot=[0,0.3])#, title='GnomView', unit='mK', format='%.2g')
hp.graticule() # print the meridians and parralels
"""


LMAX = 1024
cl = hp.anafast(m, lmax=LMAX)

P.show()


sys.exit()


