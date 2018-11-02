#!/usr/bin/python
# -*- coding: UTF-8 -*-

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from time import sleep

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

def vk2(E,x):                         # k2 = (sqrt(e-V))^2
    k2 = E-V_Q(x)
    #k2 = 6
    return k2 
			 
def numerov(n,h,k2,u):                # d2u/dt2 + k2*u(t) = 0
    b=(h**2)/12.0                         
    for i in range(1, n-1):  
       u[i+1] = (2*u[i]*(1-5*b*k2[i])-(1.+b*k2[i-1])*u[i-1])/(1+b*k2[i+1])

E  = 6 
x0 = -1
dx  = 0.01
n  = int(3*3.14/dx)
x  = [x0+i*dx for i in range(n)]
u  = [0 for i in range(n)]
k2 = [vk2(E, x[i]) for i in range(n)]

u[0] = 3.0
u[1] = 2.9999

numerov(n,dx,k2,u)

plt.plot(x, u)
plt.show()
