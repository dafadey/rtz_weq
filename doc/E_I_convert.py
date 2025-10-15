#!/usr/bin/python3

import math

e0 = 8.854 * 1e-12 #dielectric permittivity F/m

c = 3 * 1e8 # speed of light m/s

k = 1/math.sqrt(e0*c/2)

print("V/m = " +str(k)+" sqrt(W/m^2)")

print(math.sqrt(2.4*1e9/math.pi/math.pow(0.0001,2))*k)

lam = 780 * 1e-9 # 780 nm wavelength

k0 = 2 * math.pi / lam;

P_cr = 2.4 * 1e9 # 2.4GW

E_a = 5.14 * 1e11 # V/cm

K = math.pi * e0 * c * math.pow(E_a, 2) / (math.pow(k0, 2) * P_cr)

print("dimensionless Kerr koefficient is " + str(K))

# 0.01 * math.pow(0.0011744,2) = K * math.pow(E,2)

E_new = math.pow(0.01 / K, .5) * 0.0011744

print("E_new=" +str(E_new))

t0 = lam/2/math.pi/c
print('t0='+str(t0*10**15)+' fs')

me = 9.1093837139*10**(-28) # g
qe = 4.8032*10**(-10) #esu

hbar =  6.6262*10**(-27)/(2*math.pi) #erg.sec

wa = me * qe**4 / hbar**3

print('w_a='+f'{wa:e}'+' s^-1')

A_dim = t0 * wa

print('A_dim='+str(A_dim))


print('n0='+str(k0**2/4/math.pi/qe**2*me*c**2))


print('m c^2/(4*pi e^2)*k0^2='+str(me * c**2/(4*math.pi*qe**2)*k0**2))

print('xi3='+str(0.0377/E_a**2)+' [m^2/V^2]')

chi3 = 0.0377/E_a**2

print('xxx='+str(chi3*me/(4*math.pi*qe**2 * E_a)))
