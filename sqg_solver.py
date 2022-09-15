#This program solves the 2d SQG Equations
#Codes written by Ramjee Sharma
#All rights reserved
#6/11/2022
#Import Python libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import timeit
from numpy import linalg as LA #for matrix norm

#start time
start_time = timeit.default_timer()

#Define global parameters
kappa=0.01
alpha=0.3
nSteps=6000
dt=0.001

#space discritization
N=128
x0=0.0
xN=2*np.pi
x=np.linspace(x0,xN,N, endpoint=False)
y=np.linspace(x0,xN,N, endpoint=False)
x1,x2=np.meshgrid(x,y)

#Wave Numbers
k=np.fft.fftfreq(int(N))*int(N)
k1,k2=np.meshgrid(k,k)
k11=np.sqrt(k1*k1+k2*k2)
kalpha=np.power(k11,alpha)

k11[0,0]=0.001 
k=1./k11 #1/|k|
k[0,0]=999

#2/3 dealiasing
myfilter=np.ones(N)
for i in range(N):
    if i>N/3-1 and i<5*N/6-2:
        myfilter[i]=0

myfilter1,myfilter2=np.meshgrid(myfilter,myfilter)
myfilter=myfilter1*myfilter2

#Initial data
theta=np.sin(x1)*np.sin(x2)+np.cos(x2)

#L2 Norm of the initial data
print('L2 norm of initial data:',np.round(LA.norm(theta),2)) #L2 norm

#L_infinity of initial data
print('L_infinity norm of initial data:',np.round(LA.norm(theta, np.inf),2)) 

#Fourier transform
theta1=np.fft.fft2(theta)

#Gradient of Theta
def grad(theta):
    theta1=np.fft.fft2(theta)
    theta1x1=1j*k1*theta1
    theta1x2=1j*k2*theta1
    thetax1=np.real(np.fft.ifft2(theta1x1))
    thetax2=np.real(np.fft.ifft2(theta1x2))
    return np.sqrt(thetax1*thetax1+thetax2*thetax2)

#L_infinity norm of gradient theta
print('L_infinity norm of initial grad theta:',np.round(LA.norm(grad(theta), np.inf),2))

#RK4 Routine,
def rk4(x,dt):
    k1=f(x)
    k2=f(x+0.5*dt*k1)
    k3=f(x+0.5*dt*k2)
    k4=f(x+dt*k3)
    x=x+(dt/6.0)*(k1+2.0*k2+2.0*k3+k4)#Time integration
    return x

#Forward Euler Routine
def euler(x,dt):
    return x+dt*f(x)

# Right side function
def f(f1):
    f1=f1*myfilter
    u1=(-1.0)*1.0j*k2*f1*k 
    v1=1.0j*k1*f1*k 
    u=np.real(np.fft.ifft2(u1)) 
    v=np.real(np.fft.ifft2(v1))
    f=np.real(np.fft.ifft2(f1)) 
    uf=u*f
    vf=v*f
    uf1=np.fft.fft2(uf)*myfilter 
    vf1=np.fft.fft2(vf)*myfilter 
    return -(1.0j*k1*uf1+1.0j*k2*vf1+kappa*kalpha*f1)

#Time integration Loop
for i in range (nSteps):
    thetafinal=np.real(np.fft.ifft2(theta1))#theta before possible blow up
    #theta1=rk4(theta1,dt)#RK4
    theta1=euler(theta1,dt)#Forward Euler
    theta=np.real(np.fft.ifft2(theta1))
    if np.isnan(np.max(np.max(theta))):
        n=round(i*dt,6)
        print('Singularity at t = '+str(n))
        break

#Calculating CPU run time for computation
end_time = timeit.default_timer()
print('Total CPU run Time: ', np.round(end_time - start_time,2),'seconds')

#Calculating total time
n=round(nSteps*dt,4)

#Calculating L2 and L_infinity Norms
print('t= '+str(n )+ ', L2 norm of final data: '+str(np.round(LA.norm(thetafinal),2))) #L2 norm
print('t= '+str(n )+ ', L_infinity norm of final data: '+str(np.round(LA.norm(thetafinal, np.inf),2))) #Infinity Norm

#max norm of gradient theta
print('t= '+str(n )+ ', L_infinity norm of final grad theta: '+str(np.round(LA.norm(grad(thetafinal), np.inf),2)))

#Plotting contour plot of final theta
plt.contour(x1,x2,thetafinal,30)
plt.colorbar()
if kappa>0:
    plt.title('alpha = '+str(round(alpha,4)) + ',  t= '+str(n )+ ', © at Ramjee Sharma')
else:
    plt.title('t= '+str(n )+ ', © at Ramjee Sharma')
plt.show()

#Plottinng Power Spectrum 
nx = thetafinal.shape[0]
theta1 = np.fft.fftn(thetafinal)
amplitudes = np.abs(theta1)**2
k = np.fft.fftfreq(nx) * nx
k1,k2 = np.meshgrid(k, k)
knrm = np.sqrt(k1**2 + k2**2)
knrm = knrm.flatten()
amplitudes = amplitudes.flatten()
kbins = np.arange(0.5, nx//2+1, 1.)
kvals = 0.5 * (kbins[1:] + kbins[:-1])
Abins, _, _ = stats.binned_statistic(knrm, amplitudes,
                                        statistic = "mean",
                                        bins = kbins)
Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)
plt.loglog(kvals, Abins)
plt.xlabel("$k$")
plt.ylabel("$P(k)$")
if kappa>0:
    plt.title('power spectrum of SQG, alpha = '+str(round(alpha,4)) + ',  t= '+str(n )+ ', © at Ramjee Sharma')
else:
    plt.title('power spectrum of SQG, t= '+str(n )+ ', © at Ramjee Sharma')
plt.show()
