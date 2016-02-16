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
from numpy.linalg import inv
import matplotlib.pyplot as plt
import astropy
import sys

f, (P1, P2) = P.subplots(2, sharex=True, sharey=False, gridspec_kw=dict(height_ratios=[3,1]))
P1.set_title('Hubble diagram')
P1.scatter(z,mu,label='SN Ia (740)')
P1.errorbar(z,mu, linestyle='', xerr=dz,yerr=dmu,ecolor='blue',alpha=1.0,zorder=0)
P1.set_xlabel('$z$')
P1.set_ylabel('$\\mu = m_b - (M - \\alpha X_1 - \\beta C)$')
P1.set_xscale('log')
P1.set_xlim(0.008,1.2)
P1.set_ylim(30,48)
P1.plot(xfunc,yfunc)
P1.legend(bbox_to_anchor=(0.5, 0.9))
P1.text(0.1, 34, string_chi2Ndof, ha='left', rotation=0)
P1.text(0.1, 32, string_pM, ha='left', rotation=0)
P2.scatter(z,diff)
P2.errorbar(z,diff, linestyle='', xerr=dz,yerr=dmu,ecolor='blue',alpha=1.0,zorder=0)
P2.set_xlabel('$z$')
P2.set_ylabel('$\\mu - \\mu_{\\Lambda {\\rm CDM}}$')
P2.set_xscale('log')
P2.set_xlim(0.008,1.2)
P2.set_ylim(-1,1)
P2.plot(xfunc,yzero)
f.subplots_adjust(hspace=0)
P.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)




"""
Pcovar = N.array([[2.84e-4,-1.51e-5,-2.42e-4,1.75e-4,-3.51e-5],[-1.51e-5,4.91e-5,-2.35e-5,9.40e-6,-2.68e-5],[-2.42e-4,-2.35e-5,7.11e-3,3.80e-5,-7.78e-5],[1.75e-4,9.4e-6,3.8e-5,2.07e-4,-1.22e-4],[-3.51e-5,-2.68e-5,-7.78e-5,-1.23e-4,1.76e-4]])
"""
