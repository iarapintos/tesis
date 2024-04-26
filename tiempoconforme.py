'''lo que haremos ahora es considerar la variacion de masa respecto del tiempo
conforme, que se relaciona con el tiempo cosmico segun

                d eta / dt = 1/ a(eta)
                
luego, las variaciones de masa nos quedan como

  dM / d eta = dM / d eta \exp + dM / d eta \ac
  
donde las variaciones estan dadas por

    a causa de la expansion 
            
            dM/ d eta \ecp = m_b d a(eta)/d eta 
            
    a causa de la acrecion
    
            dM/ d eta \ac = (27*pi*G^2)/c^5 rho(eta) M^2(eta), 
            
            rho(eta) = 3/(8piG) 1/(a_b^2 eta_b^4) eta^2 (1+(eta/eta_b)^2)^-3

la idea es resolver las ecuaciones diferenciales por separado y luego la suma
total de las masas

'''

#%%

import matplotlib.pyplot as plt 
import numpy as np

#%% 

from astropy import constants as ast

m_sol = ast.M_sun.cgs.value   # masa del sol en gramos
G = ast.G.cgs.value           #constante de gravitacion universal en cgs
c = ast.c.cgs.value           #velocidad de la luz en cgs\

#%%

from scipy import constants as sci

pi = sci.pi

#%%

from scipy.integrate import odeint

'definimos algunas constantes varias, con a_b=1'

m_b = m_sol*10**(-7)
eta_b = 2.356e-4

'definimos la variacion de masa'
def model_bounce(M, eta):
    dM_b = (m_b/eta_b**2) * eta*(1+(eta/eta_b)**2)**(-1/2)
    return dM_b

'los puntos de tiempo conforme'
eta = np.linspace(-1,1, 1000)
eta_inicial = eta[0]

'''tomamos la condicion inicial m0, que viene de considerar que m0 y m_b se
relacionan segun 

                m0 = m_b*a(eta_inic) '''

m0 = m_b * (1+(eta_inicial/eta_b)**2)**(1/2)


'resolvemos el modelo'
sol_bounce = odeint(model_bounce, m0, eta)

# plot results
plt.plot(eta,sol_bounce)
plt.xlabel('tiempo conforme')
plt.ylabel('variacion de masa')
plt.yscale('log')
plt.show()

#%%

'graficamos como es el factor de escala dependiendo del tiempo conforme'

eta = np.linspace(-1, 1, 100)
eta_b = 2.356e-4

a = (1+(eta/eta_b)**2)**(1/2)

plt.plot(eta, a)
plt.xlabel('tiempo conforme')
plt.ylabel('factor de escala a(eta)')
plt.yscale('log')
plt.show()

#%%

'nos fijamos la relacion entre el tiempo cosmico y el tiempo conforme'
eta = np.linspace(-1, 1, 100)
eta_b = 2.356e-4

t = 1/2 * (eta_b*np.arcsinh(eta/eta_b)+eta*np.sqrt((eta/eta_b)**2+1))

plt.plot(eta, t)
plt.xlabel('tiempo conforme eta')
plt.ylabel('tiempo cosnico t')
#plt.yscale('log')
plt.show()

#%%

'verificamos como es la variacion de masa respecto a la acrecion'

eta = np.linspace(-1,1, 1000)
M = m_sol


rho = (3/8*pi*G)*(1/eta_b**4)*(eta**2)*(1+(eta/eta_b)**2)**(-3)
dM_ac = (27*pi*G**2/c**5)*rho*M**2

plt.plot(eta, dM_ac)
plt.xlabel('tiempo conforme eta')
plt.ylabel('variacion de masa debido a acrecion')
plt.yscale('log')
plt.show()

#%%

'hacemos la integracion de manera medio truca'

eta = np.linspace(-10,10, 1000)
eta_inic = eta[0]
M_inic = m_sol*(10**8)
eta_b = 2.356e-4


rho = (3/8*pi*G)*(1/eta_b**4)*(eta**2)*(1+(eta/eta_b)**2)**(-3)
dM_ac = (27*pi*G**2/c**5)*rho*M_inic**2

M = M_inic + dM_ac * (eta - eta_inic)

plt.plot(eta, M)
plt.xlabel('tiempo conforme eta')
plt.ylabel('variacion de masa debido a acrecion')
#plt.yscale('log')
plt.show()


#%%

eta = np.linspace(-1,1, 1000)
rho = (3/8*pi*G)*(1/eta_b**4)*(eta**2)*(1+(eta/eta_b)**2)**(-3)


plt.plot(eta, rho)
plt.xlabel('tiempo conforme eta')
plt.ylabel('rho(t)')
plt.yscale('log')
plt.show()

#%%

eta = 0
rho = (3/8*pi*G)*(1/eta_b**4)*(eta**2)*(1+(eta/eta_b)**2)**(-3)

print(rho)

#%%

'resolvemos ahora la variacion de masa respecto de la acrecion, sin considerar la condicion inicla de masas'

from scipy.integrate import odeint

'definimos algunas constantes varias, con a_b=1'
eta_b = 2.356e-4

'definimos la variacion de masa'
def model_ac(M, eta):
    
    rho = (3/8*pi*G)*(1/eta_b**4)*(eta**2)*(1+(eta/eta_b)**2)**(-3)

    dM_ac = (27*pi*G**2/c**5)*rho*M**2
    return dM_ac

'los puntos de tiempo conforme'
eta = np.linspace(-1,1, 1000)

'condicion inicial'
m0 = m_sol*10**2

'resolvemos el modelo'
sol_bounce = odeint(model_bounce, m0, eta)

# plot results
plt.plot(eta,sol_bounce)
plt.xlabel('tiempo conforme')
plt.ylabel('variacion de masa por acrecion')
plt.yscale('log')
plt.show()


#%%
eta_b = 2.356e-4
eta = np.linspace(-1,1, 1000)

integreal_rho = ((eta_b**2+eta**2)*np.arctan(eta/eta_b)-eta_b**3*eta+eta_b*eta**3)/8*eta_b*()


#%%

'''tomamos la condicion inicial m0, que viene de considerar que m0 y m_b se
relacionan segun 

                m0 = m_b*a(eta_inic) '''

m0 = m_b * (1+(eta_inicial/eta_b)**2)**(1/2)


'resolvemos el modelo'
sol_bounce = odeint(model_bounce, m0, eta)

# plot results
plt.plot(eta,sol_bounce)
plt.xlabel('tiempo conforme')
plt.ylabel('variacion de masa')
plt.yscale('log')
plt.show()




