from math import sqrt
from math import log10
from math import log
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


#fonction qui permet la creation d'un plot en projection Mollweide (a lire en entier)
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
    #x=-x    # reverse the scale: East to the left
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

#Boucle pour le decalage des angles (au centre)
for i in range(len(ra)):
	if ra[i]>180:
		ra[i]-=360

#using Astropy SkyCoordinates
coord=[]#coordinates in degree
skyra=[]#right asension
skydec=[]#declination
gc=[] #galactic coordinates 
gcb=[] # galactical longitude
gcl=[]#galactical latitude

coord = SkyCoord(ra=ra,dec=dec, frame='icrs', unit='deg')
skyra = coord.ra.degree
skydec = coord.dec.degree

#display of the coordinates (right asencion and declinaison) of each supernovae , the color tell which experience took the measure, legend to be done
colz=[]#colors
colz=exp

rab=[]
decb=[]
ray=[]
decy=[]
rar=[]
decr=[]
rag=[]
decg=[]

#Each supernovae, classified by experiment
for i in range(len(exp)):
	if exp[i]==3.0:
		rab.append(ra[i])
		decb.append(dec[i])
	elif exp[i]==2.0:
		rag.append(ra[i])
		decg.append(dec[i])
	elif exp[i]==1.0:
		ray.append(ra[i])
		decy.append(dec[i])
	elif exp[i]==4.0:
		rar.append(ra[i])
		decr.append(dec[i])


#Each supernovae, classified by experiment using Astropy for the Mollweide projection
blue =  SkyCoord(ra=rab,dec=decb, frame='icrs', unit='deg', )
green =  SkyCoord(ra=rag,dec=decg, frame='icrs', unit='deg', )
yellow=  SkyCoord(ra=ray,dec=decy, frame='icrs', unit='deg', )
red=  SkyCoord(ra=rar,dec=decr, frame='icrs', unit='deg', )
skyrab = blue.ra.degree
skydecb = blue.dec.degree
skyrag = green.ra.degree
skydecg = green.dec.degree
skyray = yellow.ra.degree
skydecy = yellow.dec.degree
skyrar = red.ra.degree
skydecr = red.dec.degree


#display of the coordinates (right asencion and declinaison) of each supernovae , the color tell which is the redshift associated
for i in range(len(exp)):
	if zcmb[i] >= 0.0:
		colz[i] = 'blue'
	if zcmb[i] >= 0.25:
		colz[i] = 'green'
	if zcmb[i] >= 0.5:
		colz[i] = 'yellow'
	if zcmb[i] >= 0.75:
		colz[i] = 'red'


#separation of data and color by redshifts

rabz=[]
decbz=[]
rayz=[]
decyz=[]
rarz=[]
decrz=[]
ragz=[]
decgz=[]

for i in range(len(exp)):
	#if zcmb[i] >= 0.0:
	colz[i] = 'blue'
	rabz.append(ra[i])
	decbz.append(dec[i])
	if zcmb[i] >= 0.25:
		colz[i] = 'green'
		ragz.append(ra[i])
		decgz.append(dec[i])
	if zcmb[i] >= 0.5:
		colz[i] = 'yellow'
		rayz.append(ra[i])
		decyz.append(dec[i])
	if zcmb[i] >= 0.75:
		colz[i] = 'red'
		rarz.append(ra[i])
		decrz.append(dec[i])

#display of the coordinates (right asencion and declinaison) of each supernovae , the color tell which is the redshift associated, for the Mollweide projection

blue =  SkyCoord(ra=rabz,dec=decbz, frame='icrs', unit='deg', )
green =  SkyCoord(ra=ragz,dec=decgz, frame='icrs', unit='deg', )
yellow=  SkyCoord(ra=rayz,dec=decyz, frame='icrs', unit='deg', )
red=  SkyCoord(ra=rarz,dec=decrz, frame='icrs', unit='deg', )

skyrabz = blue.ra.degree
skydecbz = blue.dec.degree
skyragz = green.ra.degree
skydecgz = green.dec.degree
skyrayz = yellow.ra.degree
skydecyz = yellow.dec.degree
skyrarz = red.ra.degree
skydecrz = red.dec.degree

#exactly the same thing for galactical coordinates
skyrabzgc = blue.galactic.l.degree
skydecbzgc = blue.galactic.b.degree
skyragzgc = green.galactic.l.degree
skydecgzgc = green.galactic.b.degree
skyrayzgc = yellow.galactic.l.degree
skydecyzgc = yellow.galactic.b.degree
skyrarzgc = red.galactic.l.degree
skydecrzgc = red.galactic.b.degree

#shift the 'ra' axis (I don't know why i should do it again...)
for i in range(len(skyrabzgc)):
	if skyrabzgc[i]>180:
		skyrabzgc[i]-=360
for i in range(len(skyragzgc)):
	if skyragzgc[i]>180:
		skyragzgc[i]-=360
for i in range(len(skyrayzgc)):
	if skyrayzgc[i]>180:
		skyrayzgc[i]-=360
for i in range(len(skyrarzgc)):
	if skyrarzgc[i]>180:
		skyrarzgc[i]-=360


#Plots

#Experiments
fig = P.figure(figsize=(10, 5))
P.scatter(rab,decb, c='blue', label = 'Low-z')
P.scatter(rag,decg, c='green', label = 'SDSS')
P.scatter(ray,decy, c='yellow', label = 'SNLS')
P.scatter(rar,decr, c='red', label = 'HST')
P.title('740 Ia supernovae classified by experiment')
P.xlabel('Right ascension in degree')
P.ylabel('Declination in degree')
P.grid()
P.legend()
P.savefig('Images plots/Map/740 Ia supernovae classified by experiment.png',dpi=180)
#P.show()
#sys.exit()

#Experiments en projection de Mollweide
fig = P.figure(figsize=(10, 5))
plot_mwd(skyrab,skydecb,color='blue', label = 'Low-z')
plot_mwd(skyrag,skydecg,color='green', label = 'SDSS')
plot_mwd(skyray,skydecy,color='yellow', label = 'SNLS')
plot_mwd(skyrar,skydecr,color='red', label = 'HST')
P.title('Angular distribution of the JLA type 1a SN in celestial coordinates, classified by experiment')
#P.xlabel('Right ascension in degree')
#P.ylabel('Declination in degree')
P.legend()
P.savefig('Images plots/Map/Angular distribution of the JLA type 1a SN in Celestial coordinates.png',dpi=180)


#Redshifts
#P.scatter(ra,dec, c=colz, label = '?')
fig = P.figure(figsize=(10, 5))
P.scatter(rabz,decbz, c='blue', label = 'z > 0.0')
P.scatter(ragz,decgz, c='green', label = 'z > 0.25')
P.scatter(rayz,decyz, c='yellow', label = 'z > 0.5')
P.scatter(rarz,decrz, c='red', label = 'z > 0.75')
P.title('ok?')
P.title('740 Ia supernovae classified by redshifts')
P.xlabel('Right ascension in degree')
P.ylabel('Declination in degree')
P.grid()
P.legend()
P.savefig('Images plots/Map/740 Ia supernovae classified by redshifts.png',dpi=180)

#Redshifts en projection de Mollweide
fig = P.figure(figsize=(10, 5))
plot_mwd(skyrabz,skydecbz,color='blue', label = 'z > 0.0')
plot_mwd(skyragz,skydecgz,color='green', label = 'z > 0.25')
plot_mwd(skyrayz,skydecyz,color='yellow', label = 'z > 0.5')
plot_mwd(skyrarz,skydecrz,color='red', label = 'z > 0.75')
P.title('Angular distribution of the JLA type 1a SN in celestial coordinates, classified by experiment')
#P.xlabel('Right ascension in degree')
#P.ylabel('Declination in degree')
P.legend()
P.savefig('Images plots/Map/Angular distribution of the JLA type 1a SN in celestial coordinates,classified by redshifts.png',dpi=180)


#Plot of each supernovae in galactical coordinates, classified by redshifts
fig = P.figure(figsize=(10, 5))
P.scatter(skyrabzgc,skydecbzgc, c='blue', label = 'z > 0.0')
P.scatter(skyragzgc,skydecgzgc, c='green', label = 'z > 0.25')
P.scatter(skyrayzgc,skydecyzgc, c='yellow', label = 'z > 0.5')
P.scatter(skyrarzgc,skydecrzgc, c='red', label = 'z > 0.75')
P.title('740 Ia supernovae in galactical coordinates, classified by redshifts')
P.xlabel('Galactical latitude in degree')
P.ylabel('Galactical longitude in degree')
P.grid()
P.legend()
P.savefig('Images plots/Map/Angular distribution of the JLA type 1a SN in ceglestial coordinates,classified by redshifts.png',dpi=180)


#Plot of each supernovae in galactical coordinates, classified by redshifts (Mollweide)
fig = P.figure(figsize=(10, 5))
plot_mwd(skyrabzgc,skydecbzgc,color='blue', label = 'z > 0.0')
plot_mwd(skyragzgc,skydecgzgc,color='green', label = 'z > 0.25')
plot_mwd(skyrayzgc,skydecyzgc,color='yellow', label = 'z > 0.5')
plot_mwd(skyrarzgc,skydecrzgc,color='red', label = 'z > 0.75')
P.title('740 Ia supernovae in galactical coordinates, classified by redshifts')
P.xlabel('Galactical latitude in degree')
P.ylabel('Galactical longitude in degree')
P.legend()
P.savefig('Images plots/Map/Angular distribution of the JLA type 1a SN in galactical coordinates,classified by redshifts.png',dpi=180)

P.show()


