#!/usr/bin/python3

import math

print('this calculates cumulative coefficietns for probablity formulae in cuda kernels from known constants')

# A(l,m,I) * (B(l,m,i)/F) ^ p(l,m,I) * exp(C(l,m,i)/F)
# F = E / E_a
# l, m, I properties of molecule
# expression in numerical algorithm:
# exp(ln(A) + (ln(B) - ln(F)) * p - C/F)
# A0 * exp((B0 - ln(F)) * p0 - C0/F)

def gamma(i) :
  if i - int(i) != 0 or i == 0 :
    print('gamma error! argument is not integer or null ' + str(i))
  res = 1
  for ii in range(1,i) :
    res *= ii
  return res

def Q(l, m) :
  f = 0
  if m <= 0 or m % 2 == 0 :
    f = 1
  else :
    f = -1
  return f * math.sqrt((2*l+1)*gamma(l+abs(m)+1)/(2*gamma(l-abs(m)+1)))

def adk(B_m, kappa, Z_c, E_t, m) :
  #pow(B_m, 2)/(pow(2, fabs(m))*tgamma(fabs(m)+1)*pow(kappa,2*Z_c/kappa-1));
  A = B_m**2 / (2**abs(m) * gamma(abs(m)+1) * kappa**(2*Z_c/kappa-1))
  #pow(2*pow(kappa,3)/F, 2*Z_c/kappa-fabs(m)-1);
  B = 2*kappa**3
  p = 2*Z_c/kappa-abs(m)-1
  #exp(-2*pow(kappa,3)/(3*F));
  C = 2/3*kappa**3
  return [A, B, p, C]

def getO2() :
  Ip=12.06 # first ionization potential for N2 molecule
  I_H=13.59 # ionization potential for hydrogen atom
  kappa = math.sqrt(Ip/I_H)
  Z_c = 1
  m = 1 # magnetic quantum number: for sigma orbitals, m=0; for pi orbitals, m=1
  # From Phys. Rev. A 81, 033423 (2010)
  C_2_m = 0.52
  C_4_m = 0.03
  # From Tong et al., Phys. Rev. A 66, 033402 (2002)
  #C_2_m = 0.62
  #C_4_m = 0.03
  B_m = C_2_m*Q(2,m)+C_4_m*Q(4,m)
	#return 1/3*w_adk(B_m, kappa, Z_c, E_t, m); # 1/3 is due to averaging over the angles, Phys. Rev. A 66, 033402 (2002)
  [A, B, p, C] = adk(B_m, kappa, Z_c , 0 , m)
  return [2*A , B, p, C]

def getN2() :
  Ip=15.58 # first ionization potential for O2 molecule
  I_H=13.59 # ionization potential for hydrogen atom
  kappa = math.sqrt(Ip/I_H)
  Z_c = 1
  m = 0 # magnetic quantum number: for sigma orbitals, m=0; for pi orbitals, m=1
  # From Phys. Rev. A 81, 033423 (2010)
  C_0_m = 2.68
  C_2_m = 1.1
  C_4_m = 0.06
  # From Tong et al., Phys. Rev. A 66, 033402 (2002)
  #C_0_m = 2.02
  #C_2_m = 0.78
  #C_4_m = 0.04
  B_m = C_0_m*Q(0,m)+C_2_m*Q(2,m)+C_4_m*Q(4,m);
  #return 2*w_adk(B_m, kappa, Z_c, E_t, m); // 2.0 is due to averaging over the angles, Phys. Rev. A 66, 033402 (2002)	
  [A, B, p, C] = adk(B_m, kappa, Z_c , 0 , m)
  return [A/3 , B, p, C]

def toKernelPrint(par) :
  print(str(par[0]*math.exp(par[2]*math.log(par[1])))+' * exp( - '+ str(par[2]) + ' * ln(F) - ' + str(par[3]) +' / F)')


#print('gamma test:')
#for i in range(0,13):
#  print(str(i+1)+'! = ' +str(gamma(i+1)))
print('-------ionization data is from-------------------')
print('V.A. Kostin, I.D. Laryushin, A.A. Silaev and N.V. Vvedenskii, "Ionization-Induced Multiwave Mixing: Terahertz Generation with Two-Color Laser Pulses of Various Frequency Ratios" PRL 117, 035003 (2016)')
print('-------------------------------------------------')
print('O2:')
toKernelPrint(getO2())
print('N2:')
toKernelPrint(getN2())

#O2: 6.018834097463142 * exp( - 0.12307858699746754 * ln(F) - 0.5573147246191912 / F)
#N2: 9.691581554890364 * exp( - 0.8679102160158225 * ln(F) - 0.8183342643785264 / F)

print('-------inozation from numerical code-------------')

def wO2(_F) :
	F=math.fabs(_F)
	[A,B,p,C] = getO2()
	return A * math.exp(p*math.log(B)) * math.exp( - p * math.log(F+1e-7) - C / (F+1e-7))
	
def wN2(_F) :
	F=math.fabs(_F)
	[A,B,p,C] = getN2()
	return A * math.exp(p*math.log(B)) * math.exp( - p * math.log(F+1e-7) - C / (F+1e-7))

tau_0 = 35 * 1e-15 # s
c = 3 * 1.e10 # cm/s
omega_0 = 2 * math.pi/(780 * 1e-7/c)
print('tau_0='+f'{tau_0:e}'+' omega_0='+f'{omega_0:e}' + ' T='+f'{2*math.pi/omega_0:e}')
A0 = 0.0733613
#A0 = 0.1344
#A0 = 0.0366807
def env(t):
	return A0 * math.exp(-t**2/2/tau_0**2)
def E(t) :
	return env(t) * math.sin(2*omega_0*math.pi*t)

def percentage(wadk) :
	n0 = (wadk(A0+.5e-7)-wadk(A0-.5e-7))/1e-7 * A0 / wadk(A0)
	wbar = (2/math.pi/n0)**.5 * wadk(A0)
	tau_i = tau_0 / n0**.5
	print("n0="+f'{n0:e}'+' tau_i='+f'{tau_i:e}')
	return (2*math.pi)**.5 * tau_i * wbar

def percentage_direct_integ(func, wadk) :
	count = 4096
	tbox = tau_0 * 6
	dt = tbox/count
	res=0	
	for i in range(0,count):
		t=-tbox*.5+dt*i
		#print('t='+f'{t}'+' E='+f'{func(t):e}')
		res += wadk(func(t))*dt
	return res

print("O2")
print('sigma_a='+f'{percentage(wO2):e}')
print('sigma_n='+f'{percentage_direct_integ(E, wO2):e}')
print("N2")
print('sigma_a='+f'{percentage(wN2):e}')
print('sigma_n='+f'{percentage_direct_integ(E, wN2):e}')

def envn(arg) :
	return (.5 + .5 * math.cos(arg)) if math.fabs(arg) < math.pi else .0

def En(t) :
	#t == 2*pi/omeg0_0 -> tt = 2*pi
	tt = t * omega_0
	arg = tt / (tau_0 * omega_0) * math.pi
	#print('t='+f'{t:e}'+' tt='+f'{tt:e}'+' arg='+f'{arg:e}')
	return A0 * envn(arg) * math.sin(tt)

w0 = 4.13e16 # s^{-1}
nAir = 2.7e19 # cm^{-3}
nO2 = .2 * nAir
nN2 = .8 * nAir
n0 = 1.8355e21

print("O2")
print('sigma_n='+f'{percentage_direct_integ(En, wO2):e}')
print("N2")
print('sigma_n='+f'{percentage_direct_integ(En, wN2):e}')

overall = (percentage_direct_integ(En, wO2) * nO2 + percentage_direct_integ(En, wN2) * nN2) * w0 / n0
print('Overall:='+f'{overall:e}')

def percentage_direct_integ_2(func, wadk) :
	count = 1024
	tbox = 250
	dt = tbox/count
	res=0	
	for i in range(0,count):
		t=-tbox*.5+dt*i
		#print('t='+f'{t}'+' E='+f'{func(t):e}')
		res += wadk(func(t))*dt
	return res

def En_2(t) :
	#t == 2*pi/omeg0_0 -> tt = 2*pi
	arg = t / (tau_0 * omega_0) * math.pi
	return A0 * envn(arg) * math.sin(t)

overall = (percentage_direct_integ_2(En_2, wO2) * nO2 + percentage_direct_integ_2(En_2, wN2) * nN2) * 17.3 / n0
print('Overall:='+f'{overall:e}')
print(nAir/n0)
