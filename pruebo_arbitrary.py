import matplotlib.pyplot as plt 
import numpy as np
import time

#%%

from astropy import constants as ast


m_sol = ast.M_sun.cgs.value   # masa del sol en gramos
G = ast.G.cgs.value           #constante de gravitacion universal en cgs
c = ast.c.cgs.value           #velocidad de la luz en cgs

#%%

from scipy.constants import physical_constants

t_planck = physical_constants["Planck time"][0]

masa_planck = physical_constants["Planck mass"][0]
m_planck = masa_planck*(10**3) #masa de planck en gramos

#%%
import mpmath

t_corrida_i = time.time()

mpmath.mp.dps = 15
mpmath.mp.pretty = True

print('entrando a la definicion del modelo')

#definimos el modelo
def model(t, M):
    if M  > 10**17:
        A_M = 5.3e25
        dMdt = -A_M / M**2
    else: 
        A_M = 7.8*10**26  
        dMdt = -A_M / M**2
    return dMdt

M_0 = (6.5*10**2)*m_planck
    
tiempo = np.logspace(np.log10(10**-25), np.log10(10**6), 100)

print('entrando a la eq diferencial')

ode_model = lambda M, t: model(t, M)

Masa = mpmath.odefun(ode_model, -10e6, M_0, tol=0.01)
    
mpmath.plot(Masa, xlim=[-10e6, 10e6], ylim=None, points=200)

t_corrida_f = time.time()

t_total = t_corrida_f - t_corrida_i

print(f'el programa tardo {t_total} segundos en correr')


#%% 

import mpmath

mpmath.mp.dps = 15
mpmath.mp.pretty = True


print('defino ctes')

M_0 = (6.5*10**2)*m_planck
A_M = 7.8e26

print('entrando a la eq diferencial')

tiempo = np.logspace(-25, 6, 10)
print(tiempo)

Masa = mpmath.odefun(lambda M, t: -A_M/M**2, 2e-25, M_0, tol=0.01)

#for t in tiempo:
#   print(Masa(t))
datos = []

for t in [2e-25, 2.7e-22, 7.7e-19, 1.5e-15, 7.4e-8, 4.2e-2, 15, 1.5e2]:
    m_t = Masa(t)
    if m_t>0:
        datos.append(m_t)

print(datos)

#vemos que claramente no funciona bien con masas muy pequenas, el bh se 
#evapora muy muy rapido 



#plt.plot(tiempo, y)
#plt.xlabel('Tiempo [t]')
#plt.ylabel('Masa [M]')
#plt.title('Evolución de la masa en función del tiempo')
#plt.show()

#%%
import mpmath

mpmath.mp.dps = 15
mpmath.mp.pretty = True


M_0 = 2.5e28
A_M = 5.3e25


tiempo = np.logspace(-25, 0, 30)

Masa = mpmath.odefun(lambda M, t: -A_M/M**2, 2e-25, m_sol*10e6, tol=0.01)

y = [Masa(t) for t in tiempo]
 
print(y)
    

plt.plot(tiempo, y)
plt.xlabel('tiempo')
plt.ylabel('masa')

#plt.xlim(1e-25, 1e-9)
plt.xscale('log')
plt.yscale('log')

plt.show()

#%%

masas = np.logspace(np.log10((6.5*10**2)*m_planck), np.log10((6.5*10**15)*m_planck), 8)

t = np.logspace(np.log10(10**-25), np.log10(10**6), 25)

resultados = []

for M_0 in masas:     
    if M_0 > 10**17:
        A_M = 5.3*10**25
        M = (M_0**3 - 3*A_M*t)**(1/3)
        M_pos = [x if x > 0 else 0 for x in M]
       
    else: 
        A_M = 7.8*10**26
        M = (M_0**3 - 3*A_M*t)**(1/3)
        M_pos = [x if x > 0 else 0 for x in M]
    
    resultados.append(M_pos)
    
for M in resultados:
    plt.plot(t, M)

plt.xlabel('tiempo')
plt.ylabel('masa')

#plt.xlim(1e-25, 1e-9)
plt.xscale('log')
plt.yscale('log')

plt.show() 


#%%

import mpmath

mpmath.mp.dps = 15
mpmath.mp.pretty = True


print('defino ctes')

masas = (10**(-13)*m_sol, m_sol*10**(-5), m_sol*10**6) 

for M_0 in masas: 
    if M_0 > 1e17: 
        A_M = 5.3e25
        Masa = mpmath.odefun(lambda M, t: -A_M/M**2, 2e-25, M_0, tol=0.01)
        


M_0 = (6.5*10**2)*m_planck
A_M = 7.8e26

print('entrando a la eq diferencial')

tiempo = np.logspace(-25, 6, 10)
print(tiempo)

Masa = mpmath.odefun(lambda M, t: -A_M/M**2, 2e-25, M_0, tol=0.01)

#for t in tiempo:
#   print(Masa(t))
datos = []

for t in [2e-25, 2.7e-22, 7.7e-19, 1.5e-15, 7.4e-8, 4.2e-2, 15, 1.5e2]:
    m_t = Masa(t)
    if m_t>0:
        datos.append(m_t)

print(datos)

