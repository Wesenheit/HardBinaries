import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import math
import sys
class particle:

    def __init__(self,mass):
        self.name=self
        self.vec=np.array([0,0,0,0,0,0])
        self.mass=mass #in M_sun

    def defvec(self,x,y,z,vx,vy,vz):
        self.vec=np.array([x,y,z,vx,vy,vz])  #AU,AU,AU,km/s,km/s,km/s
    
    def addvec(self,X):
        self.vec=self.vec+X

def force(name1,name2):
    r=math.sqrt((name1.vec[0]-name2.vec[0])**2+(name1.vec[1]-name2.vec[1])**2+(name1.vec[2]-name2.vec[2])**2)
    ax=-0.0059686133264*name2.mass*(name1.vec[0]-name2.vec[0])/r**3
    ay=-0.0059686133264*name2.mass*(name1.vec[1]-name2.vec[1])/r**3 
    az=-0.0059686133264*name2.mass*(name1.vec[2]-name2.vec[2])/r**3 
    a=[ax,ay,az]
    return np.array(a)

def D(X,aw):        #derivativ of vector
    k=1/(1.495e8)   #constant changing km/s to AU/s
    W=[k*X[3],k*X[4],k*X[5],aw[0]/1000,aw[1]/1000,aw[2]/1000]
    return np.array(W)

def rungekutta(name1,dt,aw):       # fourth order runge kutta methode
    k1=D(name1.vec,aw)
    k2=D(name1.vec+dt*k1/2,aw)
    k3=D(name1.vec+dt*k2/2,aw)
    k4=D(name1.vec+dt*k3,aw)
    deltax=(dt/6)*(k1+k4+2*(k3+k2))
    return deltax
#33333333333333333333
a=2   #au
e=0.1
m1=1 #Msun
m2=1 #Msun
a1=a*m2/(m2+m1)
a2=a*m1/(m1+m2)
vsum=29.47*math.sqrt((1+e)/(1-e))*math.sqrt((m1+m2)/a)
v1=vsum*m2/(m2+m1)
v2=vsum*m1/(m2+m1)
star1=particle(m1)
star2=particle(m2)
star1.defvec(a1,0,0,0,v1,0)
star2.defvec(-a2,0,0,0,-v2,0)
T=a**(3/2)*math.sqrt(m1+m2)
###################
c=sys.argv[4]
intruder=particle(float(sys.argv[1]))
intruder.defvec(float(c),20,0,0,-1*float(sys.argv[2]),0)
Ep=1/2*intruder.mass*(intruder.vec[4]**2)


#######3
dt=3600
Tmax=365*24*3600*T*5
t=0
X1=[[],[],[]]
Y1=[[],[],[]]
Z=[[],[]]
while t<Tmax:
    delta1=rungekutta(star1,dt,force(star1,star2)+force(star1,intruder))
    delta2=rungekutta(star2,dt,force(star2,star1)+force(star2,intruder))
    delta3=rungekutta(intruder,dt,force(intruder,star1)+force(intruder,star2))
    star1.addvec(delta1)
    star2.addvec(delta2)
    intruder.addvec(delta3)
    t=t+dt
    X1[0].append(star1.vec[0])
    X1[1].append(star1.vec[1])
    X1[2].append(star1.vec[2])
    Y1[0].append(star2.vec[0])
    Y1[1].append(star2.vec[1])
    Y1[2].append(star2.vec[2])
    Z[0].append(intruder.vec[0])
    Z[1].append(intruder.vec[1])
"""
fig = plt.figure()
ax = fig.gca(projection='3d')

plt.plot(X1[0],X1[1],X1[2])
plt.plot(Y1[0],Y1[1],Y1[2])
"""
Ek=1/2*intruder.mass*(intruder.vec[3]**2+intruder.vec[4]**2+intruder.vec[5]**2)
file=open('out.txt','a')
L=[str(intruder.mass), " ", str(Ep)," ",str(Ek)," ",str(c),'\n']
file.writelines(L)
file.close()
#plt.xlim(-8,8)
#plt.ylim(-8,8)
plt.plot(X1[0],X1[1])
plt.plot(Y1[0],Y1[1])
plt.plot(Z[0],Z[1])
plt.savefig('plot{}.png'.format(str(sys.argv[3])))



