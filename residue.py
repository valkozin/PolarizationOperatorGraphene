from scipy import integrate, real, imag
import numpy as np
from math import sqrt, pi, cos, sin, exp
import cmath

def residue(f, z0, r):
    '''
    :return: residue of f(z) at z=z0 by calculating the complex integral over a circle of radius r
    '''
    integrand = lambda phi, r: f(z0 + r * cmath.exp(1j * phi)) * 1j * r * cmath.exp(1j * phi)
    integrand_Re = lambda phi, r: real(integrand(phi, r))
    integrand_Im = lambda phi, r: imag(integrand(phi, r))
    integral_Re = integrate.quad(integrand_Re, -pi, pi, args=(r,))
    integral_Im = integrate.quad(integrand_Im, -pi, pi, args=(r,))

    return integral_Re[0] + 1j * integral_Im[0] / (2 * pi * 1j), integral_Re[1:], integral_Im[1:]

r = 0.1
z0 = 0

func = lambda z: 1/z

print(residue(func, z0, r)[0])








