#%% 

#voy a probar intentar resolver la eq diferencial general con el 
#metodo anterior pero la verdad es que no le tengo nada de fe

#%% 

import matplotlib.pyplot as plt 
import numpy as np

#%% 

from astropy import constants as ast


m_sol = ast.M_sun.cgs.value   # masa del sol en gramos
G = ast.G.cgs.value           #constante de gravitacion universal en cgs
c = ast.c.cgs.value           #velocidad de la luz en cgs

#%% 
from scipy import constants as sci
from scipy.constants import physical_constants

t_planck = physical_constants["Planck time"][0]
pi = sci.pi

masa_planck = physical_constants["Planck mass"][0]
m_planck = masa_planck*(10**3) #masa de planck en gramos

#%% 

#definimos nuestras constantes 
x_b = 9*10**37                 # constante adimensional que me define a a_b con x_b < 10^38 
a_b = 1/x_b                     # constante del bounce
T_b = t_planck*10**25          # 10^3 < T_b / t_planck < 10^40 [s]    asi que elijo uno intermedio
A_M = (5.3*10**25)              # en unidades de g^3 s^-1
lamb = 1.1056*10**(-56)        # constante cosmologica en cm^-2 
rho_0 = (lamb*c**2)/(8*pi*G)   # densidad de "hoy" en cgs. consideramos la densidad de vacio, que es la que domina


#%%

import mpmath

mpmath.mp.dps = 15
mpmath.mp.pretty = True

# constantes varias
w = 1/3
aux = 1+3*w*c**2
A = (aux**(aux/(2*w*c**2)))/(4*w**(3/2)*c**3)
rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))
P = w*(c**2)*rho

#masa inicial para el caso de la dinamica del espacio tiempo
m_0 = m_sol*(1e-10) 

#condicion inicial de masa para la resolucion de la edo
M_0 = m_sol*(1e-15) 

#las ecuaciones diferenciales 
dot_M_b = m_0*((2*a_b*t)/(3*(1-w)*T_b))*(1+(t/T_b)**2)**((-2-3*w)/(3*(1-w))) 
dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P)

tiempo = np.logspace(-25, 6, 10)

Masa = mpmath.odefun(lambda M, t: dot_M_b+dot_M_ac, 2e-25, m_sol*10e6, tol=0.01)

y = [Masa(t) for t in tiempo]

plt.plot(tiempo, y)
plt.xlabel('tiempo')
plt.ylabel('masa')

#plt.xlim(1e-25, 1e-9)
plt.xscale('log')
plt.yscale('log')

plt.show()


