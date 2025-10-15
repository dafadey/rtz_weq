import math

def r(R0, k0, z):
    return math.sqrt(R0*R0 + 4*z*z/(k0*k0*R0*R0))

def rp(R0, k0, z):
    return R0 * math.sqrt(R0*R0 + 4 * z*z/(k0*k0*R0*R0)) * math.sqrt(k0/(2 * z))


def _r(R0, z):
    return math.sqrt(R0*R0 + 64*z*z/(R0*R0))

def _rp(R0, z):
    return R0 * math.sqrt(R0*R0 + 64 * z*z/(R0*R0)) / math.sqrt(8 * z)


lmbd = 780. * math.pow(10,-9) *100 #cm 
k0 = 2.*math.pi/lmbd

R0 = 0.02 * math.pow(10,-3) #mm

x0 = 1/(math.sqrt(2)*k0)

z0 = 2/k0

beamDia = 5. * 1.7 #mm

focalLength = 85. #mm

d0 = .1 * beamDia / x0
d02 = d0 * d0

z_foc = .1 * focalLength / z0

dw = math.sqrt(d02/2. - math.sqrt(d02*d02/4.-1024.*z_foc*z_foc))

z_shift=13402/1.7 #.4/z0*160./dw

aaa = dw/2
r_0 = math.sqrt(math.pow(z_shift/aaa,2.0)*64.0 + math.pow(aaa,2.0))

print("z_shift="+str(z_shift)+" : " +str(z_shift*z0)+"cm")
print("r_0 = "+str(r_0)+" : "+str(r_0*x0)+"cm")
print("dw = "+str(dw)+" : " + str(dw*x0)+"cm")
print("aaa = "+str(aaa)+" : " + str(aaa*x0)+"cm")
print("d0 = " + str(d0) + " : "+str(d0*x0)+"cm")

#print(r_0*2/d0*0.015)

print(math.pow(_r(aaa,z_foc)/aaa*0.0011,2)*0.01)
