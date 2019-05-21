from scipy import integrate, real, imag
import numpy as np
from math import sqrt, pi, cos, sin, exp
import cmath
print(integrate.quad(lambda x: sin(x), -100, 100, weight='cauchy', wvar=0))
quit()
sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]])

v = 1
k_max = np.inf

def G(omega, kx, ky, s, delta):
    return 1/(2*sqrt(kx**2+ky**2)*(omega-s*v*sqrt(kx**2+ky**2)+1j*s*delta))*(sqrt(kx**2+ky**2)*np.identity(2)+s*(kx*sx+ky*sy))

def ff1(k1x, k1y, s1, omega2, k2x, k2y, s2):
    return 1/(2*sqrt(k1x**2+k1y**2))*(sqrt(k1x**2+ky**2)*np.identity(2)+s1*(k1x*sx+k1y*sy))*\
           1/(2*sqrt(k2x**2+k2y**2)*(omega2-s2*v*sqrt(k2x**2+k2y**2)))*(sqrt(k2x**2+k2y**2)*np.identity(2)+s2*(k2x*sx+k2y*sy))

def ff2(omega1, k1x, k1y, s1, k2x, k2y, s2):
    return 1/(2*sqrt(k1x**2+k1y**2)*(omega1-s1*v*sqrt(k1x**2+k1y**2)))*(sqrt(k1x**2+k1y**2)*np.identity(2)+s1*(k1x*sx+k1y*sy))*\
           1/(2*sqrt(k2x**2+k2y**2))*(sqrt(k2x**2+k2y**2)*np.identity(2)+s2*(k2x*sx+k2y*sy))

def lim_phi1(k1, omega, kx, ky, eta):
    return [0, 2*pi]

def lim_k1(omega, kx, ky, eta):
    return [0, k_max]

def opts_phi1(k1, omega, kx, ky, eta):
    return {}

def opts_k1(omega, kx, ky, eta):
    return {}

def integrand(phi1, k1, omega, kx, ky, eta):
    #res1 = -1j*1/(2*pi)**3*k1*np.trace(ff1(k1*cos(phi1), k1*sin(phi1), -1, omega+(-v*k1), kx+k1*cos(phi1), ky+k1*sin(phi1), 1))
    #res2 = -1j*1/(2*pi)**3*k1*np.trace(ff2(-omega-v*sqrt((kx+k1*cos(phi1))**2+(ky+k1*sin(phi1))**2), k1*cos(phi1), k1*sin(phi1), 1, kx+k1*cos(phi1), ky+k1*sin(phi1), -1))
    #return 4*(2*pi*1j*res1+2*pi*1j*res2)
    return 4*(-1j)*1/(2*pi)**3*k1**(-v)*(k1+sqrt((k1*cos(phi1)+kx)**2+(k1*sin(phi1)+ky)**2))/(v**2*(k1+sqrt((k1*cos(phi1)+kx)**2+(k1*sin(phi1)+ky)**2))**2-omega**2)
integrand_Re = lambda phi1, k1, omega, kx, ky, eta: real(integrand(phi1, k1, omega, kx, ky, eta))
integrand_Im = lambda phi1, k1, omega, kx, ky, eta: imag(integrand(phi1, k1, omega, kx, ky, eta))

omega = 0
kx = 1
ky = 0
eta = 0.0
#print(integrand(2*pi, k1, omega, kx, ky, eta))
#quit()
integral_Re = integrate.nquad(integrand_Re, [lim_phi1, lim_k1], args=(omega, kx, ky, eta), opts=[opts_phi1, opts_k1])
integral_Im = integrate.nquad(integrand_Im, [lim_phi1, lim_k1], args=(omega, kx, ky, eta), opts=[opts_phi1, opts_k1])
print(integral_Re[0]+1j*integral_Im[0], integral_Re[1:], integral_Im[1:])
