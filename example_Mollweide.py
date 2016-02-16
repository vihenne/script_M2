from math import *
import pylab as P
import numpy as N 
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
import sys


#Ouverture du fichier de donnees des 740 supernovae:
datafile=open('/data/software/jla_likelihood_v6/data/jla_lcparams.txt','r')

#Lecture de la 1ere ligne (noms des colonnes) puis passage a la 2eme:
dataline=datafile.readline()
dataline=datafile.readline()
#   Formation d'une liste composee des elements de la 2eme ligne (1ere supernova) ; "split" separe les donnees des differentes colonnes de la ligne :
liste=dataline.split()
#   Listes a representer graphiquement & initialisation des Variables
zcmb=[]
mB=[]
dmB=[]
X1=[]
dX1=[]
C=[]
dC=[]
M_stell=[]
ra=[]
dec=[]
exp=[]

c=299792.458
H=0.000070
omgM=0.295
alpha=0.141
beta=3.101
M_b=-19.05
delta_M=-0.070

#Boucle -> remplissage des tableaux de donnees extraites des colonnes du tableau des 740 supernovae
while dataline != "" :
#"append" permet d'ajouter un element a la liste a laquelle il est applique (ex: zmcb, mB...), "float" permet de convertir une chaine de nombre (telle que lue dans le tableau des 740 supernovae) et de la convertir en float. 
	zcmb.append(float(liste[1]))
	mB.append((float(liste[4])))
	dmB.append(float(liste[5]))
	X1.append(float(liste[6]))
	dX1.append(float(liste[7]))
	C.append(float(liste[8]))
	dC.append(float(liste[9]))
	M_stell.append(float(liste[10]))
	exp.append(float(liste[17]))
	ra.append(float(liste[18]))
	dec.append(float(liste[19]))
#passage a la ligne suivante du tableau :       
	dataline=datafile.readline()
	liste=dataline.split()


#fonction qui permet la creation d'un plot en projection Mollweide
def plot_mwd(RA,Dec,org=0,title='Mollweide projection', projection='mollweide',color='black',label='label'):
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
    ax = fig.add_subplot(111, projection=projection,)
    ax.scatter(N.radians(x),N.radians(Dec),c=color, label = label)  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)


coord = SkyCoord(ra=ra,dec=dec, frame='icrs', unit='deg')
skyra = coord.ra.degree
skydec = coord.dec.degree


#RA - DEC
fig = P.figure(figsize=(10, 5))
P.scatter(ra,dec)
P.title('740 Ia supernovae')
P.xlabel('Right ascension in degree')
P.ylabel('Declination in degree')
P.grid()
P.legend()


#Mollweide
fig = P.figure(figsize=(10, 5))
plot_mwd(skyra,skydec)
P.title('Angular distribution of the JLA type 1a SN in celestial coordinates')
P.xlabel('Right ascension in degree')
P.ylabel('Declination in degree')
P.legend()
P.show()

