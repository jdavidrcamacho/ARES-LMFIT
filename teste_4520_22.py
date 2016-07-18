# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 12:36:24 2016

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

data= np.loadtxt('tmp_4520_22.txt',dtype=None)
y=data[:,1]
x=data[:,0]
#figure()
#pl.plot(data[:,0],data[:,1])
#pl.plot(x,y)

#ver onde fazemos  corte
#LIMITS SPECTRA FOR FIT: lli: 4537.560 -- llf: 4538.260  <-do ARES
a0=np.where(x<4519.551)
a1=list(a0[-1])
a2=a1[-1]
b0=np.where(x<4520.691)
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

params = [-0.009, 4519.64, 0.05, \
            -0.005, 4519.81, 0.08, \
            -0.05, 4520.00, 0.05, \
            -0.08, 4520.24, 0.05, \
            -0.01, 4520.53, 0.05,\
            -0.003, 4520.64, 0.05]
constrains = [(-0.01,-0.006), (4519.60,4519.65), (0.06,0.04), \
                (-0.008,-0.001), (4519.80,4519.85), (0.1,0.07), \
                (-0.06,-0.04), (4519.95,4520.05), (0.08,0.02), \
                (-0.12,-0.08), (4520.20,4520.25), (0.08,0.02), \
                (-0.015,-0.009), (4520.50,4520.55), (0.08,0.02), \
                (-0.005,-0.001), (4520.60,4520.65), (0.08,0.02)]
                
# constrains do codigo = boundiries              
mod, out, init = lmfit_ngauss_constrains(xfinal,yfinal, params, constrains)  
#figure()
pl.plot(xfinal, yfinal)
pl.plot(xfinal, init+1, 'k--')
print(out.fit_report(min_correl=0.9))
pl.plot(xfinal, out.best_fit+1, 'r-')
pl.show()