'''aca calcularemos las soluciones numericas a la ecuacion diferencial de la 
densidad de agujeros negros 
'''

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

pi = sci.pi

#%%

'definimos algunas constantes varias'

eta_0 =  1.43e7  #el tiempo conforme hoy
eta_b = 2.356e-4    #escala del bounce

a_b = 1/((1+(eta_0/eta_b)**2)**(1/2))

m_0 = m_sol  #la masa del bh en el momento del bounce

#%%

'''queremos calcular cuanto vale la constante de normalizacion A

            A = rho_BH (t(eta_cc)) / (M_cc c^2)
            
entonces tenemos que encontrar rho_BH 

            rho_BH = beta' * rho_rad 

con 

            beta' = beta / (1-beta)

'''

#%%

'verificamos que el beta que estamos tomando es el correcto (cumple restricciones)'

#todos estos datos fueron sacados del paper de Josan (2009)

#usamos el beta de las mediciones del WMAP (M_BH > 5e14 g)

#condicion inicial
eta_ci = -2.6e4 #son aproximadamente los 47k anios donde domina la rad
a_ci = a_b*(1+(eta_ci/eta_b)**2)**(1/2)
M_ci = m_0 * a_ci

f_M = (1/3)**(3/2)

beta_dato = 1.5e-19 * (M_ci/(f_M * 5e14))**(1/2)

print(np.log10(beta_dato))
print(M_ci)
print(np.log10(M_ci))

# ESTE VALOR NO ESTA PERMITIDO POR LAS RESTRICCIONES

#%%

'verificamos que el beta que estamos tomando es el correcto (cumple restricciones)'

#metemos mano para ver que pasa con otro beta

#condicion inicial
eta_ci = - 2.6e4 #son aproximadamente los 300000 anios donde domina la rad
a_ci = a_b*(1+(eta_ci/eta_b)**2)**(1/2)
M_ci = m_0 * a_ci

beta_dato = 1e-15   #sacado deLgrafico de las restricciones

print(np.log10(beta_dato))
print(M_ci)
print(np.log10(M_ci))


#%%

beta_calc = beta_dato / (1-beta_dato)

#definimos la densidad de radiacion para la condicion inicial
rho_rad_ci = 3/ (8*pi*G)*eta_b**(-4) * a_b**(-2) *eta_ci**2*(1+(eta_ci/eta_b)**2)**(-3)


#la condicion inicial sera 
rho_BH_ci = beta_calc * rho_rad_ci
print(rho_BH_ci)

#%%
print(M_ci)
print(m_0)  #preguntar

#%% 

#la cte de normalizacion nos queda
A = rho_BH_ci / (M_ci * c**2)

print(A)


#%%

from scipy.integrate import odeint

'definimos la variacion de densidad de bhs respecto al tiempo conforme'

cte = A*m_0 * c**2

def densidad_bh(rho, eta):
    a = (1+(eta/eta_b)**2)**(1/2) #el factor de escala
    drho = - a_b**2 *(eta/eta_b) * (1/a) * (cte+3*(1/a)*rho)
    return drho

print(cte)

#%%

'los puntos de tiempo conforme'
eta = np.linspace(-2.6e4,2.6e4, 1000)

a = (1+(eta/eta_b)**2)**(1/2) #el factor de escala
drho = - a_b**2 *(eta/eta_b) * (1/a) * (cte+3*(1/a)*rho_BH_ci)

plt.plot(eta,drho)
plt.xlabel('eta [s]')
plt.ylabel('rho_bh/d eta [erg/cm^3/s]')
plt.yscale('log')
plt.xscale('symlog')
plt.show()

#%%
'los puntos de tiempo conforme'
eta = np.linspace(-2.6e4,2.6e4, 1000)

'resolvemos el modelo'
sol_densidad_bh = odeint(densidad_bh, rho_BH_ci, eta)

# plot results
plt.plot(eta,sol_densidad_bh)
plt.xlabel('eta [s]')
plt.ylabel('rho_bh [erg/cm^3]')
plt.yscale('symlog')
plt.xscale('symlog')
plt.show()


#%% 

'suponete que la densidad es muuy despreciable'

cte = A*m_0 * c**2

def densidad_bh(rho, eta):
    a = (1+(eta/eta_b)**2)**(1/2) #el factor de escala
    drho = - a_b**2 *(eta/eta_b) * (1/a)* cte
    return drho

'los puntos de tiempo conforme'
eta = np.linspace(-2.6e4,2.6e4, 1000)

'resolvemos el modelo'
sol_densidad_bh = odeint(densidad_bh, rho_BH_ci, eta)

# plot results
plt.plot(eta,sol_densidad_bh)
plt.xlabel('eta [s]')
plt.ylabel('rho_bh [erg/cm^3]')
plt.yscale('symlog')
plt.xscale('symlog')
plt.show()


#%% 

'suponete que la constante es muuy despreciable'

cte = A*m_0 * c**2

def densidad_bh(rho, eta):
    a = (1+(eta/eta_b)**2)**(1/2) #el factor de escala
    drho = - 3* a_b**2 *(eta/eta_b) * (1/a**2) * rho
    return drho

'los puntos de tiempo conforme'
eta = np.linspace(-2.6e4,2.6e4, 1000)

'resolvemos el modelo'
sol_densidad_bh = odeint(densidad_bh, rho_BH_ci, eta)

# plot results
plt.plot(eta,sol_densidad_bh)
plt.xlabel('eta [s]')
plt.ylabel('rho_bh [erg/cm^3]')
plt.yscale('symlog')
plt.xscale('symlog')
plt.show()

#%%

from scipy.integrate import solve_ivp

cte = A*m_0 * c**2

def densidad_bh(rho, eta):
    a = (1+(eta/eta_b)**2)**(1/2) #el factor de escala
    drho = - a_b**2 *(eta/eta_b) * (1/a) * (cte+3*(1/a)*rho)
    return drho

'los puntos de tiempo conforme'
eta = np.linspace(-2.6e4,2.6e4, 1000)
#eta = np.linspace(-1e3,1e3, 1000)

'resolvemos el modelo'
sol_densidad_bh_ivp = solve_ivp(densidad_bh, (eta[0], eta[-1]), [rho_BH_ci], t_eval=eta, method='LSODA')

dens_BH = sol_densidad_bh_ivp.y[0]

# plot results
plt.plot(eta,sol_densidad_bh_ivp.y[0])
plt.xlabel('eta [s]')
plt.ylabel('rho_bh [erg/cm^3]')
#plt.yscale('symlog')
plt.xscale('symlog')
plt.savefig('densidad_BH_preliminar.png')
plt.show()

#%%

'ahora queremos sacar la densidad del fluido cosmologico'


'los puntos de tiempo conforme'
eta = np.linspace(-2.6e4,2.6e4, 1000)
#eta = np.linspace(-1e3,1e3, 1000)


cte_2 = (3*c**2*a_b**4)/(8*pi*G)

a = (1+(eta/eta_b)**2)**(1/2) #el factor de escala

rho_CF = cte_2 *(eta**2/eta_b**2)*(1/a**6)-dens_BH

# plot results
plt.plot(eta, rho_CF)
plt.xlabel('eta [s]')
plt.ylabel('rho_cf [erg/cm^3]')
#plt.yscale('symlog')
plt.xscale('symlog')
plt.show()

#%%

'ahora calculamos el parametro de la eq de estado para el fluido cosmologivo'

eta = np.linspace(-2.6e4,2.6e4, 1000)
a = (1+(eta/eta_b)**2)**(1/2) #el factor de escala

'defunimos la derivada de rho_CF respecto del tiempo conforme'
cte_3 = (3*c**2*a_b**4)/(8*pi*G*eta_b**2)
cte = A*m_0 * c**2
drho_BH = - a_b**2 *(eta/eta_b) * (1/a) * (cte+3*(1/a)*dens_BH)

drho_CF = cte_3 * (2*eta * (1/a**6)-(6*a_b**2/eta_b)*eta**3*(1/a**8)) - drho_BH

'definimos Q'
Q = -cte * a_b**2 *(eta/eta_b)* (1/a**2)

'el parametro nos queda'
omega_CF = -1 -(1/3)*(eta_b/a_b**2) * (a**3/eta) * (1/rho_CF) * (drho_CF*(1/a)+Q)


# plot results
plt.plot(eta, omega_CF)
plt.xlabel('eta [s]')
plt.ylabel('omega_CF []')
#plt.yscale('symlog')
plt.xscale('symlog')
plt.show()












