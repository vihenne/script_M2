from __future__ import division
from scipy.integrate import quad
from math import sqrt
from math import log10
from math import log
import pylab as P
import numpy as N 
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches
from scipy.optimize import fmin
from operator import add

#   Ouverture du fichier de donnees des 740 supernovae:
datafile=open('/Users/rosnet/Documents/LSST/JLA2014/jla_likelihood_v6/data/jla_lcparams.txt','r')


#   Lecture de la 1ere ligne (noms des colonnes) puis passage a la 2eme:
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
c=299792.458
H=0.000070
omgM=0.295
alpha=0.141
beta=3.101
M_b=-19.05
delta_M=-0.070


#Fonction donnant les mu theoriques sous forme d'un tableau (fonction de omega_matter et la distance-luminosite)
def intfun(y,x):
    return 1/sqrt(((1+y*x)*((1+x)**2))-(1-y)*x*(2+x))
def fitfundL(z,omgM):
    MU=[]
    for i in range (len(z)): 
        zz=z[i]
        MU.append(dL_z(zz,omgM))
    return MU  
def dL_z(zz,omgM): 
    mu_zz = 5*log10((1+zz)*c*(quad(intfun,0,zz,args=(omgM))[0]/(10*H)))
    return mu_zz  

    
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
#passage a la ligne suivante du tableau :       
      dataline=datafile.readline()
      liste=dataline.split()

#Fonction de construction des muexp sous forme d'un tableau (formule du papier de 2014)
def muexp(mB,X1,C,alpha,beta,M_b,delta_M,M_stell):
    mu=[]
    for i in range(len(mB)):
        if M_stell[i]<10:
            mu.append(mB[i]-M_b+alpha*X1[i]-beta*C[i])
        else : 
            mu.append(mB[i]-M_b-delta_M+alpha*X1[i]-beta*C[i])
    return mu

#Fonction de construction des dmuexp sous forme d'un tableau (formule de propagation des incertitudes)
def dmuexp(dmB,dX1,dC,alpha,beta):
    dmu=[]
    for i in range(len(dmB)):
        dmu.append(sqrt(dmB[i]**2+(alpha*dX1[i])**2+(beta*dC[i])**2))
    return dmu

#Construction des mu et dmu a l'aide des fonctions muexp,dmuexp
mu=muexp(mB,X1,C,alpha,beta,M_b,delta_M,M_stell)
dmu=dmuexp(dmB,dX1,dC,alpha,beta)
        

#Fonction de Fit "Curve_fit" a un parametre Omega_matter ; elle retourne la valeur du parametre optimale (popt) ainsi que la matrice de covariance (pcov)
print '*** Fit with 1 free parameter'
omG = [0.295]
popt, pcov = curve_fit(fitfundL,zcmb,mu,omG,dmu,True)
perr = N.sqrt(N.diag(pcov))
print 'Omega_M=', popt[0], ' +-',perr[0]

xfunc = N.linspace(0.001,2,1000)
yfunc = N.zeros(len(xfunc))
omgM = popt[0]
#Representation de la courbe theorique avec le parametre omgm obtenu :
yfunc = fitfundL(xfunc,omgM)


#Fonction de construction du Chi2 (avec p=[omgm,alpha,beta,M_b,delta_M]= les 5 parametres a faire fitter)
def chi2(p,z,mB,dmB,X1,dX1,C,dC,M_stell):
    omgM = p[0]
    alpha = p[1]
    beta = p[2]
    M_b = p[3]
    delta_M = p[4]
    chi2tot=0
    for i in range (len(z)):
        zz = z[i]
        if M_stell[i]<10:
            mu_z = mB[i]-M_b+alpha*X1[i]-beta*C[i]
        else : 
            mu_z = mB[i]-M_b-delta_M+alpha*X1[i]-beta*C[i]
        dmu_z = sqrt(dmB[i]**2+(alpha*dX1[i])**2+(beta*dC[i])**2)
        chi2tot = chi2tot + ((mu_z-dL_z(zz,omgM))/dmu_z)**2
    return chi2tot
#calcul du chi2 initial
chi2curfit=chi2([omgM,alpha,beta,M_b,delta_M],zcmb,mB,dmB,X1,dX1,C,dC,M_stell)
print 'chi2 =',chi2curfit

#Minimisation du chi2 grace a "fmin", en fonction des 5 parametres de fit (retournant les parametres ideaux, mais sans les incertitudes associees)
print '*** Fit with 5 free parameters'
p0=[omgM,alpha,beta,M_b,delta_M]
p = fmin(chi2,p0,args=(zcmb,mB,dmB,X1,dX1,C,dC,M_stell))

#Fonction de construction de la matrice de covariance des 5 parametres (en vue d'obtenir les incertitudes sur les 5 parametres)
from numpy.linalg import inv
def covar(param,z,mB,dmB,X1,dX1,C,dC,M_stell):
# "covar" est initialisee comme matrice de taille 5x5 (nb de parametres x nb de parametres) remplie de zeros
  covar=N.zeros(shape=(len(param),len(param)))
  P1=N.zeros(shape=(len(param)))
  P2=N.zeros(shape=(len(param)))
  P3=N.zeros(shape=(len(param)))
  E=N.zeros(shape=(len(param)))
#l'erreur est fixee a 1% de la valeur des parametres ; cette valeur peut etre controlee et modifiee
  err=0.01
  for i in range(len(param)):
    P1[i]=param[i]
    P2[i]=param[i]
    P3[i]=param[i]
  for i in range(len(param)):
    E[i]=param[i]*err
    P1[i]=P1[i]*(1+err)
    for j in range(len(param)):
      if i==j:
        P2[j]=P2[j]*(1+2*err)
#calcul de la derivee partielle double du chi2 par rapport au parametre i=j
        covar[i][j]=((chi2(P2,z,mB,dmB,X1,dX1,C,dC,M_stell)-2*chi2(P1,z,mB,dmB,X1,dX1,C,dC,M_stell)+chi2(param,z,mB,dmB,X1,dX1,C,dC,M_stell))/(2*E[i]*E[j]))
        P2[j]=P2[j]/(1+2*err)
      else:      
        P1[j]=P1[j]*(1+err)
        P2[i]=P2[i]*(1+err)
        P3[j]=P3[j]*(1+err)
        E[j]=param[j]*err
#calcul de la derivee partielle croisee du chi2 par rapport aux parametres i et j
        covar[i][j]=((chi2(P1,z,mB,dmB,X1,dX1,C,dC,M_stell)-chi2(P2,z,mB,dmB,X1,dX1,C,dC,M_stell)-chi2(P3,z,mB,dmB,X1,dX1,C,dC,M_stell)+chi2(param,z,mB,dmB,X1,dX1,C,dC,M_stell))/(2*E[i]*E[j]))
        P1[j]=P1[j]/(1+err)
        P2[i]=P2[i]/(1+err)
        P3[j]=P3[j]/(1+err)
    P1[i]=P1[i]/(1+err)
#la matrice construite est l'inverse de la matrice de covariance, il convient donc de l'inverser : 
  covar=inv(covar)
  return covar

#Calcul de la matrice de covariance 
Pcovar= covar(p,zcmb,mB,dmB,X1,dX1,C,dC,M_stell)
print Pcovar
#Calcul des incertitudes sur les 5 parametres grace a la racine carree des #valeurs de la diagonale de la matrice de covariance 
Errors=N.zeros(5)
Errors=N.sqrt(N.diag(Pcovar))
print Errors


#Affichage des valeurs de fit des 5 parametres et de leurs incertitudes
print 'Omega_M=', p[0], '+-', Errors[0]
print 'alpha=', p[1], '+-', Errors[1]
print 'beta=', p[2], '+-', Errors[2]
print 'M_b=', p[3], '+-', Errors[3]
print 'delta_M=', p[4], '+-', Errors[4]

#Nouveau calcul des valeurs des mu et des dmu avec les parametres ajustes
mufinal=muexp(mB,X1,C,p[1],p[2],p[3],p[4],M_stell)
dmufinal=dmuexp(dmB,dX1,dC,p[1],p[2])


#Fonction de calcul des ecarts entre mu_exp et mu_theorique sous forme de tableaux
def ecarts(z,mu,omgM):
    ecarts=[]
    for i in range(len(z)):
        zz = z[i]
        ecarts.append(mu[i]-dL_z(zz,omgM))
    return ecarts

#Calcul de ces ecarts entre experience et theorie
ecarts1=ecarts(zcmb,mu,omgM)
ecarts5=ecarts(zcmb,mufinal,p[0])
        

#Plots

#Premier scatter de points effectue sans corrections des effets (C, M_b, X1) et sans incertitudes, simplement mb en fonction de zcmb
P.figure()
P.scatter(zcmb, mB,s=3)
P.xlabel('zcmb')
P.ylabel('mb')
P.xscale('log')
P.title('Premier scatter')
#Premier fit a 1 parametre suite a l'utilisation de "Curve_fit"
P.figure()
P.scatter(zcmb,mu,s=1)
P.errorbar(zcmb,mu, linestyle='', xerr=0,yerr=dmu,ecolor='blue',alpha=1.0,zorder=0)
P.xlabel('zcmb')
P.ylabel('mu')
P.xscale('log')
P.plot(xfunc,yfunc,'r-')
P.title('Premier fit a 1 parametre')
#Fit a 5 parametres apres minimisation du chi2
P.figure()
P.plot(xfunc,fitfundL(xfunc,p[0]),'r-')
P.errorbar(zcmb,mufinal,linestyle='', xerr=0,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
P.xlabel('zcmb')
P.ylabel('mu')
P.xscale('log')
P.title('Fit a 5 parametres')
#Ecarts entre mu_exp et mu_theorique dans le fit a 1 parametre
P.figure()
P.scatter(zcmb,ecarts1,s=1)
P.errorbar(zcmb,ecarts1, linestyle='', xerr=0,yerr=dmu,ecolor='blue',alpha=1.0,zorder=0)
P.xscale('log')
P.xlabel('zcmb')
P.ylabel('muexp-muth')
P.title('Ecarts - fit a 1 parametre')
#Ecarts entre mu_exp et mu_theorique dans le fit a 5 parametres
P.figure()
P.scatter(zcmb,ecarts5,s=1)
P.errorbar(zcmb,ecarts5, linestyle='', xerr=0,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
P.xscale('log')
P.xlabel('zcmb')
P.ylabel('muexp-muth')
P.title('Ecarts - fit a 5 parametres')

P.show()

