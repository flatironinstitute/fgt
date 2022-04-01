from numpy import *
from pylab import *

aa = loadtxt('fort.33')



x = aa[:,0]
y = aa[:,1]
u = aa[:,2]
v = aa[:,3]

u[abs(u)<1e-16] = 1e-16
v[abs(v)<1e-16] = 1e-16
erru = abs(u-v)
erru[abs(erru)<1e-16] = 1e-16
u = log(abs(u))/log(10.0)
v = log(abs(v))/log(10.0)
erru = log(abs(erru))/log(10.0)

figure(1)
clf()
tricontourf(x,y,u)
colorbar()

figure(2)
clf()
tricontourf(x,y,v)
colorbar()

figure(3)
clf()
tricontourf(x,y,erru)
colorbar()


bb = loadtxt('fort.34')
xt = bb[:,0]
yt = bb[:,1] 

ut = bb[:,2]
vt = bb[:,3]



ut[abs(ut)<1e-16] = 1e-16
vt[abs(vt)<1e-16] = 1e-16
errut = abs(ut-vt)
errut[abs(errut)<1e-16] = 1e-16
ut = log(abs(ut))/log(10.0)
vt = log(abs(vt))/log(10.0)
errut = log(abs(errut))/log(10.0)

figure(4)
clf()
tricontourf(xt,yt,ut)
colorbar()

figure(5)
clf()
tricontourf(xt,yt,vt)
colorbar()


figure(6)
clf()
tricontourf(xt,yt,errut)
colorbar()
show()
