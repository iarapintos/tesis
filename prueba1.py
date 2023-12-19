import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy import constants as sci
from scipy.constants import physical_constants
from astropy import constants as ast

#ahora queremos solamente la contribucion de la acrecion 

#en principio quiero probar que pasa si no le pongo la condicion de masas de flor

#ctes grales
masa_planck = physical_constants["Planck mass"][0]
m_planck = masa_planck*(10**3) #masa de planck en gramos
m_sol = ast.M_sun.cgs.value   # masa del sol en gramos
G = ast.G.cgs.value           #constante de gravitacion universal en cgs
c = ast.c.cgs.value           #velocidad de la luz en cgs
pi = sci.pi


#constantes varias 
lamb = 1.1056*10**(-56)        # constante cosmologica en cm^-2 
rho_0 = (lamb*c**2)/(8*pi*G)   #densidad de HOY que NO es la correcta
a_b = 1
T_b = 1


#definimos el modelo 
def model_acretion(t, M): 
    w = 1/3
    aux = 1+3*w*c**2
    A = (aux**(aux/(2*w*c**2)))/(4*w**(3/2)*c**3)
    rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))
    P = w*(c**2)*rho
    
    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P)
    
    return dot_M_ac


#la escala de tiempo
t = np.linspace(-1e4, 1e4, 1000)
t_inic = t[0]

#condicion inicial de la masa (SIN RELACION)
y0 = m_sol*1e4

#resolvemos la edo
sol_acretion = solve_ivp(model_acretion, [t[0], t[-1]], [y0], method='RK45', t_eval=t)
    
#las masas 
y_acretion = sol_acretion.y[0]

#fialmente ploteamos
plt.plot(t, y_acretion, label='acrecion')
plt.xlabel('tiempo')
plt.ylabel('masa')
plt.legend()

plt.yscale('log')
#lt.xscale('log')
plt.show() 










#%% 

# con condicion de masa de flor

#ctes grales
masa_planck = physical_constants["Planck mass"][0]
m_planck = masa_planck*(10**3) #masa de planck en gramos
m_sol = ast.M_sun.cgs.value   # masa del sol en gramos
G = ast.G.cgs.value           #constante de gravitacion universal en cgs
c = ast.c.cgs.value           #velocidad de la luz en cgs
pi = sci.pi


#constantes varias 
lamb = 1.1056*10**(-56)        # constante cosmologica en cm^-2 
rho_0 = (lamb*c**2)/(8*pi*G)   #densidad de HOY que NO es la correcta
a_b = 1
T_b = 1
w = 1/3

#definimos el modelo 
def model_acretion(t, M): 
    aux = 1+3*w*c**2
    A = (aux**(aux/(2*w*c**2)))/(4*w**(3/2)*c**3)
    rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))
    P = w*(c**2)*rho
    print('imprimimos la densididad')
    print(rho)
    print('imprimimos la presion')
    print(P)


    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P)
    
    print('imprimimos la variacion')
    print(dot_M_ac)
    
    return dot_M_ac


#la escala de tiempo
t = np.linspace(-1e4, 1e4, 50)
t_inic = t[0]

#condicion inicial de la masa 
m_b = 2.5e20
y0 =  m_b *(1+t_inic**2)**((1/3)*(1-w))

#resolvemos la edo
sol_acretion = solve_ivp(model_acretion, [t[0], t[-1]], [y0], method='RK45', t_eval=t)
    
#las masas 
y_acretion = sol_acretion.y[0]

#fialmente ploteamos
plt.plot(t, y_acretion, label='acrecion')
plt.xlabel('tiempo')
plt.ylabel('masa')
plt.legend()

plt.yscale('log')
#lt.xscale('log')
plt.show() 

#%%

#numeros lindos

#ctes grales

G = 6.67         #constante de gravitacion universal en cgs
c = 300          #velocidad de la luz en cgs
pi = sci.pi


#constantes varias 
lamb = 1.      # constante cosmologica en cm^-2 
rho_0 = (lamb*c**2)/(8*pi*G)   #densidad de HOY que NO es la correcta
a_b = 1
T_b = 1
w = 1/3

#definimos el modelo 
def model_acretion(t, M): 
    aux = 1+3*w*c**2
    A = (aux**(aux/(2*w*c**2)))/(4*w**(3/2)*c**3)
    rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))
    P = w*(c**2)*rho
    #print('imprimimos la densididad')
   # print(rho)
   # print('imprimimos la presion')
  #  print(P)


    dot_M_ac = 4*pi*A*(G**2)*(c**(-5))*M**2*(c**2*rho-P)
    
  #  print('imprimimos la variacion')
   # print(dot_M_ac)
    
    return dot_M_ac


#la escala de tiempo
t = np.linspace(-1e4, 1e4, 50)
t_inic = t[0]

#condicion inicial de la masa 
m_b = 2
y0 =  m_b *(1+t_inic**2)**((1/3)*(1-w))

#resolvemos la edo
sol_acretion = solve_ivp(model_acretion, [t[0], t[-1]], [y0], method='RK45', t_eval=t)
    
#las masas 
y_acretion = sol_acretion.y[0]
print(y_acretion)

#%%

#fialmente ploteamos
plt.plot(t, y_acretion, label='acrecion')
plt.xlabel('tiempo')
plt.ylabel('masa')
plt.legend()

plt.yscale('log')
#lt.xscale('log')
plt.show() 

#%% 

# viendo que forma tiene la edo con numeros lindos 

#ctes grales
G = 6.67           #constante de gravitacion universal en cgs
c = 300        #velocidad de la luz en cgs
pi = sci.pi


#constantes varias 
lamb = 1        # constante cosmologica en cm^-2 
rho_0 = (lamb*c**2)/(8*pi*G)   #densidad de HOY que NO es la correcta
a_b = 1
T_b = 1
w = 1/3

#definimos el modelo 
def model_acretion(t, M): 
    aux = 1+3*w*c**2
    A = (aux**(aux/(2*w*c**2)))/(4*w**(3/2)*c**3)
    rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))
    P = w*(c**2)*rho
    
    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P)
    
    return dot_M_ac


#la escala de tiempo
t = np.linspace(-1e4, 1e4, 1000)
t_inic = t[0]

#condicion inicial de la masa 
m_b = 2.5e20
y0 =  m_b *(1+t_inic**2)**((1/3)*(1-w))

#resolvemos la edo
sol_acretion = solve_ivp(model_acretion, [t[0], t[-1]], [y0], method='RK45', t_eval=t)
    
if sol_acretion.success:
    y_acretion = sol_acretion.y[0]
    print(y_acretion)
    plt.plot(t, y_acretion, label='acrecion')
    plt.xlabel('tiempo')
    plt.ylabel('masa')
    plt.legend()
    plt.yscale('log')
    plt.show()
    
else: 
    print("Error en la integración.")
    print(sol_acretion.message)


#%% 


t = np.linspace(-1e4, 1e4, 40)

M = 2.5e4

lamb = 1.1056*10**(-56)        # constante cosmologica en cm^-2 
rho_0 = (lamb*c**2)/(8*pi*G)   #densidad de HOY que NO es la correcta
a_b = 1
T_b = 1

w = 1/3
aux = 1+3*w*c**2 
A = (aux**(aux/(2*w*c**2)))/(4*w**(3/2)*c**3) 
rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))  
P = w*(c**2)*rho  
dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P)

print('imprimimos la densididad')
print(rho)
print('imprimimos la presion')
print(P)
print('imprimimos la variacion')
print(dot_M_ac)






