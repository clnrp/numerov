#!/usr/bin/python
# -*- coding: UTF-8 -*-

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from time import sleep
import math

# d2R/dr2 + k2*R = 0
# k2 = c*[E-V(x)]

def V_Q(x):                                  
    if(x  > 3 and x < 5): v = 7.1                          
    else:                 v = 0
    return v

def V(x):                                  
    if (abs(x)<=500):	v = -0.001                          
    else:               v = 0
    return v

def vk2(E,x):                        # k2 = (sqrt(e-V))^2
    #k2 = E-V_Q(x)
    k2 = 6
    return k2 

def vk2_aton(x,n,l):                 # para u(r) = r*R(r) e rho = r/(a0*n)     
    rho0 = 2*n
    k2 = -(1 - rho0/x + l*(l+1)/x**2)
    return k2 
			 
def numerov(n,h,k2,u):               # d2u/dt2 + k2*u(t) = 0
    b=(h**2)/12.0                         
    for i in range(1, n-1):  
       u[i+1] = (2*u[i]*(1-5*b*k2[i])-(1.+b*k2[i-1])*u[i-1])/(1+b*k2[i+1])

def euler(n,h,k2,u,du0):             # d2u/dt2 = -k2*u(t)
    U=u[0]
    dU=du0
    for i in range(1,n-1):
        d2U=-k2[i-1]*u[i-1]   #U''
        dU=dU+d2U*h           #U'
        U=U+dU*h              #U
        u[i]=U

def runge_kutta4(n,h,k2,u,du0):
    U=u[0]
    dU=du0
    for i in range(1,n-1):
        k1=h*2

def atom():
    n = 3
    l = 2
    a0 = 0.53e-10
    x0 = 0.0001
    dx  = 0.01
    np  = int(10/dx)
    x  = [x0+i*dx for i in range(np)] # rho = r/(a0*n)
    u  = [0 for i in range(np)]
    k2 = [vk2_aton(x[i],n,l) for i in range(np)]

    u[0] = x[0]**(l+1) # u(p) ~ p^(l+1) quando p -> 0
    u[1] = x[1]**(l+1)

    for i in range(np): # r/a0 = rho*n
        x[i]=x[i]*n

    #para n = 1     
    #R0 = 2    
    #dRdx = -1 
    
    #para n = 2
    #R0 = 0.707 para n=2
    #dRdx = -0.707
    #R1 = R0+dx*dRdx
    #u[0] = x[0]*R0 # u(0) = r*R(0)
    #u[1] = x[1]*R1

    numerov(np,dx,k2,u)
    
    # calcular R = u(r)/r
    R=[u[i]/x[i] for i in range(np)]

    # densidade de probabilide
    p=[R[i]**2 * x[i]**2 for i in range(np)]
    
    # normalização
    norm=0
    for i in range(np):
        norm += R[i]**2 * x[i]**2 * dx
    
    for i in range(np):
        R[i] = R[i]/math.sqrt(norm)

    norm=0
    for i in range(np):
        norm += R[i]**2 * x[i]**2 * dx
    print norm

    # valor esperada
    r_esp=0
    for i in range(np):
        r_esp += R[i]**2 * x[i]**3 * dx
    print r_esp

    plt.plot(x, R)
    plt.plot(x, p)

    #du0  = (u[1]-u[0])/dx
    #euler(n,dx,k2,u,du0)
    #u2=u
    #R=[u2[i]/x[i] for i in range(n)]
    #plt.plot(x, R)

    plt.show()
    return R

def pPot():
    x0 = -1
    dx  = 0.01
    n  = int(3*3.14/dx)
    x  = [x0+i*dx for i in range(n)]
    u  = [0 for i in range(n)]
    #k2 = [vk2(E, x[i]) for i in range(n)]

    u[0] = 3.0
    u[1] = 2.9999

    numerov(n,dx,k2,u)
    plt.plot(x, u)

    #du0  = (u[1]-u[0])/dx
    #euler(n,dx,k2,u,du0)
    #u2=u
    #plt.plot(x, u2)
    plt.show()

R=atom()
