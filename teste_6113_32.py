# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 14:03:04 2016

@author: camacho
"""

import numpy as np
from lmfit.models import GaussianModel, ExponentialModel
import sys
import matplotlib.pyplot as plt
from numpy import sqrt, pi, exp, linspace, loadtxt
from astropy.io import fits

#------------------------------------------------------------------------------
def gaussian(x, amp, cen, wid):
  "1-d gaussian: gaussian(x, amp, cen, wid)"
  return (amp/(sqrt(2*pi)*wid)) * exp(-(x-cen)**2 /(2*wid**2))
#------------------------------------------------------------------------------

def lmfit_ngauss(x,y, *params):
  params = params[0]
  mods = []
  prefixes = []
  for i in range(0, len(params), 3):
    pref = "g%02i_" % (i/3)
    gauss_i = GaussianModel(prefix=pref)

    if i == 0:
      pars = gauss_i.guess(y, x=x)
    else:
      pars.update(gauss_i.make_params())

    A = params[i]
    l_cen = params[i+1]
    sigma = params[i+2]

    pars[pref+'amplitude'].set(A)
    pars[pref+'center'].set(l_cen)
    pars[pref+'sigma'].set(sigma)

    mods.append(gauss_i)
    prefixes.append(pref)

  mod = mods[0]

  if len(mods) > 1:
    for m in mods[1:]:
      mod += m

  print mod

  init = mod.eval(pars, x=x)
  out = mod.fit(y, pars, x=x)
  return mod, out, init
#------------------------------------------------------------------------------

def lmfit_ngauss_constrains(x,y, params, constrains):
  #params = params[0]
  #constrains = constrains[0]
  mods = []
  prefixes = []
  for i in range(0, len(params), 3):
    pref = "g%02i_" % (i/3)
    gauss_i = GaussianModel(prefix=pref)

    if i == 0:
      pars = gauss_i.guess(y, x=x)
    else:
      pars.update(gauss_i.make_params())

    A = params[i]
    limA = constrains[i]
    l_cen = params[i+1]
    limL = constrains[i+1]
    sigma = params[i+2]
    limS = constrains[i+2]

    pars[pref+'amplitude'].set(A, min=limA[0], max=limA[1])
    pars[pref+'center'].set(l_cen, min=limL[0], max=limL[1])
    pars[pref+'sigma'].set(sigma, min=limS[0], max=limS[1])

    mods.append(gauss_i)
    prefixes.append(pref)

  mod = mods[0]

  if len(mods) > 1:
    for m in mods[1:]:
      mod += m

  init = mod.eval(pars, x=x)
  out = mod.fit(y, pars, x=x)
  return mod, out, init
#------------------------------------------------------------------------------
import pylab as pl

data= np.loadtxt('tmp_6113_32.txt',dtype=None)
y=data[:,1] #inverti o sinal
x=data[:,0]
#figure()
#pl.plot(data[:,0],data[:,1])
#pl.plot(x,y)

#ver onde fazemos  corte
a0=np.where(x<6112.749)
a1=list(a0[-1])
a2=a1[-1]
b0=np.where(x<6113.429)
b1=list(b0[-1])
b2=b1[-1]
#corte
xfinal=x ;yfinal=y
for i in range(0,a2+1):
    yfinal[i]=1
for j in range(b2+1,len(y)):
    yfinal[j]=1       
#figure()
#pl.plot(xfinal,yfinal)  

### TENTATIVA DE FIT
#gaussian(x,amplitude,centro,width/largura)

#x = xfinal
#y = gaussian(xfinal, -0.1, 4534.0, 0.2) + gaussian(xfinal, -0.1, 4534.2, 0.2)

params = [-0.015, 6112.93, 0.10, \
        -0.001, 6113.13, 0.10, \
        -0.015, 6113.33, 0.08]
constrains = [(-0.02,-0.01), (6112.90,6112.95), (0.01,0.2), \
            (-0.002,-0.0005), (6113.10,6113.15), (0.01,0.02), \
            (-0.02,-0.01), (6113.30,6113.35), (0.05,0.1)]
#yfinal=yfinal+1                
# constrains do codigo = boundiries              
mod, out, init = lmfit_ngauss_constrains(xfinal,yfinal, params, constrains)  
figure()
pl.plot(xfinal, yfinal)
pl.plot(xfinal, init+1, 'k--')
print(out.fit_report(min_correl=0.9))
pl.plot(xfinal, out.best_fit+1, 'r-')
pl.show()