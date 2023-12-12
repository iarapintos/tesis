'''el verdadero pasado en limpio. aca junto todo lo que hicimos hasta ahora
y comento el por que de todo lo que hicimos


el objetivo de todas estas cuentas es calcular la variacion de masa del agujero 
negro para calcular su area y consecuentemente calcular lo que conocemos como 
filling factor, que es un factor que nos habla del tamanio del agujero negro 
en comparacion con el espacio-tiempo

'''

#%%


# importamos todos los paquetes necesarios que usaremos a lo largo del programa 

import csv
import matplotlib.pyplot as plt 
import numpy as np
import os
import pandas as pd
from scipy.integrate import solve_ivp

#%% 

from astropy import constants as ast
from astropy import units as u

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

# definimos la funcion que nos escribe el csv 
def guardar_datos_variacion(t, dot_M_b, dot_M_ac, dot_M_rh, archivo):
    datos_var = archivo + '_variacion_masa.csv'
    archivo_existe = False
    try:
        with open(datos_var, 'r') as archivo_csv:
            archivo_existe = True
    except FileNotFoundError: 
        archivo_existe = False
        
    modo_apertura = 'a' if archivo_existe else 'w'
    with open(datos_var, mode=modo_apertura, newline='') as archivo_csv:
        writer = csv.writer(archivo_csv)
        if not archivo_existe:
            writer.writerow(['tiempo', 'dot_M_b', 'dot_M_ac', 'dot_M_rh'])
        writer.writerow([t, dot_M_b, dot_M_ac, dot_M_rh])
        
#%%
        
masas = (10**(-13)*m_sol, m_sol*10**(-5), m_sol*10**6) #masas en g
#el limite superior lo saque de la tesis de edu. el resto los elegi yo para
#tener una nocion de como es el comportamiento de la variacion de estas masas

for i in range(len(masas)): 
    for j in range(len(masas)):
        tiempo = np.logspace(-25, 6, 1000, endpoint=False) 
        M = masas[i]
        m_0 = masas[j]
        nombre = str(i)+str(j)
        
        if M >= m_0:
            for t in tiempo:
                if t < -T_b or t > T_b: 
                    w = 1/3  #radiacion 
                    aux = 1+3*w*c.value**2
                    A = (aux**(aux/(2*w*c.value**2)))/(4*w**(3/2)*c.value**3)
                    rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))
                    P = w*(c**2)*rho
        
                    dot_M_b = m_0*((2*a_b*t)/(3*(1-w)*T_b))*(1+(t/T_b)**2)**((-2-3*w)/(3*(1-w))) 
                    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P) 
                    if M > 10**17*u.g:
                        A_M = (5.3*10**25)        # en unidades de g^3 s^-1
                        dot_M_rh = - A_M / M**2
                    else: 
                        A_M = (7.8*10**26)       # en unidades de g^3 s^-1
                        dot_M_rh = - A_M / M**2
                else:
                    w = 1/3   
                    aux = 1+3*w*c.value**2
                    A = (aux**(aux/(2*w*c.value**2)))/(4*w**(3/2)*c.value**3)
                    rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))  
                    P = w*(c**2)*rho                             
        
                    dot_M_b = m_0*((2*a_b*t)/(3*(1-w)*T_b))*(1+(t/T_b)**2)**((-2-3*w)/(3*(1-w))) 
                    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P) 
                    if M > 10**17*u.g:
                        A_M = (5.3*10**25)       # en unidades de g^3 s^-1
                        dot_M_rh = - A_M / M**2
                    else: 
                        A_M = (7.8*10**26)       # en unidades de g^3 s^-1
                        dot_M_rh = - A_M / M**2
        
            
                datos = guardar_datos_variacion(t, dot_M_b.value, dot_M_ac.value, dot_M_rh.value, nombre)
                
#%%

#en esta parte del codigo nos dedicamos a graficar con los datos que guardamos antes           

ruta = "/home/iara/faq/tesis"
archivos = os.listdir(ruta)

datos_csv = [archivo for archivo in archivos if archivo.endswith('.csv')]

for archivo in datos_csv: 
    datos = pd.read_csv(archivo)

    i = int(archivo[0])
    j = int(archivo[1])

    M = masas[i]
    m_0 = masas[j]
    
    t = datos['tiempo']
    dot_M_b = datos['dot_M_b']
    dot_M_ac = datos['dot_M_ac']
    dot_M_rh = -datos['dot_M_rh'] #el menos va porque si no no me deja pasarlo a escala log porque es negativo
    
    plt.xlabel('tiempo')
    plt.ylabel('variacion')
    plt.yscale('log')  
    plt.xscale('log') 
 
    plt.plot(t, dot_M_b, label = "Dinámica espacio-tiempo con m_0 = {:.1e}".format(m_0.value), color = 'turquoise')
    plt.plot(t, dot_M_ac, label = "Acrecion", color='darkorange')
    plt.plot(t, dot_M_rh , label = "Radiacion", color = 'deeppink')
    plt.axvline(x = T_b, color = 'mediumpurple', label = 'escala del bounce')
    plt.legend()
    plt.title('Variacion considerando una masa M= {:.1e}'.format(M.value))
    plt.show()
    
#%% 

''' 
ahora que ya tenemos las contribuciones, nos ponemos a resolver estas variaciones
de masa. 

primero resolvemos las ecuaciones diferenciales por separado. esto nos da una idea
de como se tiene que comportar la masa. ademas, al juntar todas las contribuciones
veremos a cual se parece. 

empezaremos a calcular la mas sencilla, que es la variacion de masa por parte
de radiacion de hawking

'''

#defino el modelo, con las correspondientes constantes
    
def model_hawking(t, M): 
    if M> 1e17:     
        A_M = 5.3e25  #en unidades de g^3 s^-1
        dMdt_h = -A_M / M**2
    else: 
        A_M = 7.8e26  #en unidades de g^3 s^-1
        dMdt_h = -A_M / M**2
    return dMdt_h
    
#definimos m_b, que es la masa que tiene el agujero negro en el momento
#del bounce. esto lo hacemos porque nos modifica la condicion inicial 
m_b = 2.5e20

# el tiempo en que integraremos
t = np.linspace(-1e6, 1e6, 1000)
t_inic = t[0]

# condicion inicial de las masas
w = 1/3
y0 =  m_b *(1+t_inic**2)**((1/3)*(1-w))
    
#resolvemos la eq diferencial
sol_haw = solve_ivp(model_hawking, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

#aqui tenemos la variacion de masa    
y_haw = sol_haw.y[0]

#finalmente ploteamos
plt.plot(t,y_haw, label= 'rad hawking') 
plt.xlabel('tiempo')
plt.ylabel('masa')
plt.legend()

#plt.yscale('log')
#plt.xscale('log')
plt.show()

    
#%% 

''' ahora resolvemos la variacion de masa a causa del espcio tiempo

 '''
    
#definimos la masa que tendra el agujero negro en el momento del bounce
#y el parametro de nuestro fluido
m_b = 2.5e20
w = 1/3
    
#definimos el modelo 
def model_bounce(t,M): 
    a_b = 1
    T_b = 1
   
    num = a_b*2*(1-w)* t * (1+(t/T_b)**2)**((-2/3)-(w/3))
    den = 3*T_b**2
   
    dot_M_b = m_b*(num/den) 
    
    return dot_M_b
    

#la escala de tiempo
t = np.linspace(-1e6, 1e6, 1000)
t_inic = t[0]

# condicion inicial de las masas, que vemos que depende de m_b
y0 =  m_b *(1+t_inic**2)**((1/3)*(1-w))
    
#resolvemos la edo
sol_bounce = solve_ivp(model_bounce, [t[0], t[-1]], [y0], method='RK45', t_eval=t)
    
#las masas 
y_bounce = sol_bounce.y[0]

#fialmente ploteamos
plt.plot(t, y_bounce, label='dinamica del espacio-tiempo')
plt.xlabel('tiempo')
plt.ylabel('masa')
plt.legend()

#plt.yscale('log')
#plt.xscale('log')
plt.show() 
    
#%% 

''' ahora queremos calcular la variacion de masa por efecto de la acrecion
del agujero negro 


'''
    
    
#%% 

''' para probar como es la mezcla de dos contribuciones, resolvemos ambas al 
mismo tiempo. para esto, hacemos 

'''
    
m_b = 2.5e10
w = 1/3

def model_hawking(t, M): 
    if M> 1e17:     
        A_M = 5.3e25     #en unidades de g^3 s^-1
        dMdt_h = -A_M / M**2
    else: 
        A_M = 7.8e26     #en unidades de g^3 s^-1
        dMdt_h = -A_M / M**2
    return dMdt_h

def model_bounce(t,M): 
    a_b = 1
    T_b = 1
   
    num = a_b*2*(1-w)* t * (1+(t/T_b)**2)**((-2/3)-(w/3))
    den = 3*T_b**2
   
    dot_M_b = m_b*(num/den) 
    
    return dot_M_b


# aqui consideramos la eq diferencial dM/dt = dot_M_b + dot_M_rh
def model_gral(t, M):
    
    a_b = 1
    T_b = 1
   
    num = a_b*2*(1-w)* t * (1+(t/T_b)**2)**((-2/3)-(w/3))
    den = 3*T_b**2
    
    dot_M_b = m_b*(num/den)
    
    if M> 1e17:     
        A_M = 5.3e25
        dot_M_rh = -A_M / M**2
    else: 
        A_M = 7.8e26
        dot_M_rh = -A_M / M**2

    
    dMdt = dot_M_b + dot_M_rh
    return dMdt

# escala de tiempo
t = np.linspace(-1e6, 1e6, 1000)
t_inic = t[0]

# condicion inicial de las masas 
y0 =  m_b *(1+t_inic**2)**((1/3)*(1-w))



# resolvemos las edos
sol_gral = solve_ivp(model_gral, [t[0], t[-1]], [y0], method='RK45', t_eval=t)
sol_haw = solve_ivp(model_hawking, [t[0], t[-1]], [y0], method='RK45', t_eval=t)
sol_bounce = solve_ivp(model_bounce, [t[0], t[-1]], [y0], method='RK45', t_eval=t)

#y las variaciones de masa 
y_gral = sol_gral.y[0] 
y_haw = sol_haw.y[0]
y_bounce = sol_bounce.y[0]

# ploteamos
plt.plot(t,y_gral, label= 'solucion general')
plt.plot(t,y_haw, label= 'rad hawking')
plt.plot(t, y_bounce, label='dinamica del espacio-tiempo')
plt.xlabel('tiempo')
plt.ylabel('masa')
plt.legend()
plt.title('Variacion considerando una masa inicial M_0= {:.1e}'.format(y0))

#plt.yscale('log')
#plt.xscale('log')
plt.show()
    
#%% 

''' ahora probamos con otras contribuciones 


'''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    