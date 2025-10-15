import math
lam=780*math.pow(10,-9)*100 #cm^-1
k0=2*math.pi/lam
z0=2/k0
x0=1/math.sqrt(2)/k0
alpha=1.8962 #gaussian beam factor
n2=4*math.pow(10,-23)*math.pow(100,2) #cm^2/W
Pcr=alpha*math.pow(lam,2)/4/math.pi/n2
print(str(Pcr/math.pow(10,9))+" GW")

print(0.0005/x0)
print(str(2231.91*x0*10)+" mm")