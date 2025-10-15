#!/usr/bin/python3
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy
from scipy.fft import fft, ifft, fftfreq
import math

#functions are from ioniz_probability_precalc

w0 = 4.13e16 # s^{-1}

def getO2(_F):
  F=math.fabs(_F)+1e-13
  return w0*6.018834097463142 * math.exp( - 0.12307858699746754 * math.log(F) - 0.5573147246191912 / F)

def getN2(_F):
  F=math.fabs(_F)+1e-13
  return w0*9.691581554890364 * math.exp( - 0.8679102160158225 * math.log(F) - 0.8183342643785264 / F)


#pulse from numerical simulation:
tau_0 = 70 * 1e-15 # s
c = 3 * 1.e10 # cm/s
Ea_Vcm = 5.14e9 # V/cm
Ea_Vm = Ea_Vcm*1e2 # V/m
omega_0 = 2 * math.pi/(780 * 1e-7/c)
A0max = .256765
A0min = .0427941
Kerr = 0.0377
t0_AU = 2.417e-17
n_crit = 1.83e21 # cm^-3

def envn(arg) :
  return (.5 + .5 * math.cos(arg)) if math.fabs(arg) < math.pi else .0

def En(A, t) :
  #t == 2*pi/omeg0_0 -> tt = 2*pi
  tt = t * omega_0
  arg = tt / (tau_0 * omega_0) * math.pi
  #print('t='+f'{t:e}'+' tt='+f'{tt:e}'+' arg='+f'{arg:e}')
  return A * envn(arg) * math.sin(tt)

def En3(A, t) :
  #t == 2*pi/omeg0_0 -> tt = 2*pi
  tt = t * omega_0 * 3
  arg = tt / (tau_0 * omega_0) * math.pi
  #print('t='+f'{t:e}'+' tt='+f'{tt:e}'+' arg='+f'{arg:e}')
  return A * envn(arg) * math.sin(tt)


nmax = 2.23 * 1e19 #cm^-3
nO2max = .2*nmax
nN2max = .8*nmax


#integrate plasma density
def n(A) :
  count = 2048*2
  tbox = tau_0 * 3
  dt = tbox/count
  res=0
  wnO2=0
  wnN2=0
  for i in range(0,count):
    t=-tbox*.5+dt*i
    wnO2 += getO2(En(A,t)) * dt
    wnN2 += getN2(En(A,t)) * dt
  nO2 = nO2max*(1-math.exp(-wnO2))
  nN2 = nN2max*(1-math.exp(-wnN2))
  return nO2 + nN2

#filter
def filter3w(omega_0, tbox, count, src):
  dw = 2*math.pi/tbox
  wi = omega_0*3/dw
  dwi_2 = omega_0/dw/2*1.3
  _srcimg = fft(src)
  srcimg = _srcimg.copy()
  ws = numpy.zeros(count//8)
  for i in range(0,count//8):
    ws[i] = i * dw
  
  for i in range(1,count//2):
    f = math.exp(-((i-wi)/dwi_2)**8)
    srcimg[i] *= f
    srcimg[count-i] *= f
  
  #fig, ax = plt.subplots()
  #ax.plot(ws, 1e-3+numpy.abs(_srcimg[0:count//8]), color='red')
  #ax.plot(ws, 1e-3+numpy.abs(srcimg[0:count//8]), color='green')
  #ax.set_yscale('log')
  #plt.show()
  return ifft(srcimg).real

def filter1w(omega_0, tbox, count, src):
  dw = 2*math.pi/tbox
  wi = omega_0/dw
  dwi_2 = omega_0/dw/2*1.3
  _srcimg = fft(src)
  srcimg = _srcimg.copy()
  ws = numpy.zeros(count//8)
  for i in range(0,count//8):
    ws[i] = i * dw
  
  for i in range(1,count//2):
    f = math.exp(-((i-wi)/dwi_2)**8)
    srcimg[i] *= f
    srcimg[count-i] *= f
  
  return ifft(srcimg).real

def energy(tbox, count, field) :
  dt = tbox/count
  nrg=0
  for i in range(0,count) :
    nrg += field[i]**2 * dt
  return nrg


# get plasma source
# n * E
def plasmaSource(A, dbg = False) :
  count = 2048*2
  tbox = tau_0 * 3
  dt = tbox/count
  wnO2=0
  wnN2=0
  n = numpy.zeros(count)
  E = numpy.zeros(count)
  ts = numpy.zeros(count)
  for i in range(0,count):
    ts[i] = i * dt
    t=-tbox*.5+dt*i
    _e = En(A,t)
    wnO2 += getO2(_e) * dt
    wnN2 += getN2(_e) * dt
    nO2 = nO2max*(1-math.exp(-wnO2))
    nN2 = nN2max*(1-math.exp(-wnN2))
    n[i] = nO2+nN2
    E[i] = _e
  src = (n * E) / nmax
  res = filter3w(omega_0, tbox, count, src) 
  if dbg == True :
    fig, ax = plt.subplots()
    ax.plot(ts, n/10, color='black')
    ax.plot(ts, E*n[count-1], color='gray')
    ax.plot(ts, src * nmax, color='red')
    ax.plot(ts, res * nmax, color='green')
    plt.show()
  return energy(tbox, count, res) / t0_AU

#print(plasmaSource(.111265, True))
#print(plasmaSource(.07275, True))

def wwwSource(A, dbg = False) :
  count = 2048
  tbox = tau_0 * 3
  dt = tbox/count
  E = numpy.zeros(count)
  ts = numpy.zeros(count)
  for i in range(0,count):
    ts[i] = i * dt
    t=-tbox*.5+dt*i
    E[i] = En(A,t)
  #_src = (.41e-15)**2 * n_crit * Kerr * (E * E * E) /  nmax / t0_AU**.5
  _src = 8.723e-47 * (E * Ea_Vm)**3 /  nmax
  src = numpy.zeros(count)
  for i in range(1,count-1):
    src[i] = (2*_src[i] - _src[i-1] - _src[i+1]) / dt**2
  res = filter3w(omega_0, tbox, count, src)
  resw = filter1w(omega_0, tbox, count, src) 
  if dbg == True :
    fig, ax = plt.subplots()
    #ax.plot(ts, n/10, color='black')
    #ax.plot(ts, E*n[count-1], color='gray')
    #ax.plot(ts, src, color='red')
    ax.plot(ts, res, color='green')
    ax.plot(ts, resw, color='blue')
    plt.show()
  return energy(tbox, count, res) / t0_AU

#print(wwwSource(1, True))

def test() :
  count = 2048
  tbox = tau_0 * 3
  dt = tbox/count
  E = numpy.zeros(count)
  E3 = numpy.zeros(count)
  ts = numpy.zeros(count)
  for i in range(0,count):
    ts[i] = i * dt
    t=-tbox*.5+dt*i
    E[i] = En(1,t)
  
  E3 = E**3
  
  fig, ax = plt.subplots()
  ax.plot(ts, E3, color='green')
  ax.plot(ts, filter1w(omega_0, tbox, count, E3), color='red')
  ax.plot(ts, filter3w(omega_0, tbox, count, E3), color='blue')
  #ax.plot(ts, E3, color='blue')
  plt.show()

#test()

npts=33

Is = numpy.zeros(npts)
ns = numpy.zeros(npts)
sps = numpy.zeros(npts)
swwws = numpy.zeros(npts)

printed1 = False
printed2 = False
for i in range(0,npts) :
  amp = A0min + (A0max-A0min)*i/(npts-1)
  #amp to I
  I = (amp/27.49*Ea_Vcm)**2
  n_final = n(amp)
  sp = plasmaSource(amp)
  swww = wwwSource(amp)
  #print(f'{I:e}'+':'+f'{n(amp):e}')
  Is[i] = I
  ns[i] = n_final
  sps[i] = sp 
  swwws[i] = swww
  if I>1e14 and printed1 == False:
    printed1 = True
    print(f'{I:e}'+' : '+f'{sp:e}')
  if I>2e14 and printed2 == False:
    printed2 = True
    print(f'{I:e}'+' : '+f'{sp:e}')

plt.rcParams['figure.figsize'] = [10,7]

fig, [ax1, ax2] = plt.subplots(1,2)
ax1.plot(Is, ns)

ax1.set_yscale('log')
ax1.set_xscale('log')

ax2.plot(Is, sps, color='purple')
ax2.plot(Is, swwws, color='green')

ax2.set_yscale('log')
ax2.set_xscale('log')

locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
ax2.yaxis.set_major_locator(locmaj)

locmin = matplotlib.ticker.LogLocator(base=10.0, subs=numpy.arange(2, 10) * .1, numticks=100)
ax2.yaxis.set_minor_locator(locmin)


plt.show()

print(math.log(2.7))
