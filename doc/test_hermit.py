import numpy as np
import math

def A(i):
    return 2+1/(i+.5)

def C(i):
    if i==0 :
        return -4
    else:
        return -1-1/(i+.5)

def B(i):
    return -1

N=13

OPP=np.zeros((N,N))


for i in range(0,N):
    for j in range(0,N):
        if i==j :
            OPP[i][j]=A(i)
        elif j==i+1 :
            OPP[i][j]=C(i)
        elif j==i-1 :
            OPP[i][j]=B(i)

#print(OPP)

print("----------------------------")

def test(v1, v2) :
    res = .0
    for i in range(0,N):
        if i==0 :
            r = .3/.8
        else :
            r = i+.5
        res = res + r*v1[i]*v2[i]
    return res

def apply(O, v):
    res = np.zeros(N)
    for i in range(0,N):
        for j in range(0,N):
            res[i] = res[i] + O[i][j]*v[j]
    return res

diff=0

for j in range(0,3):
    for i in range(0,N):
        v1=np.zeros(N)
        v2=np.zeros(N)
        v1[i]=1
        if i+j<N:
          v2[i+j]=1
        v1_OPPv2 = test(v1, apply(OPP,v2))
        OPPv1_v2 = test(apply(OPP,v1),v2)
        print(str(v1_OPPv2) + " vs " + str(OPPv1_v2))
        diff = math.fabs(v1_OPPv2 - OPPv1_v2)
    print("---------------diff="+str(diff))
