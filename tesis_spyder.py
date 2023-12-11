# importamos todos los paquetes necesarios que usaremos a lo largo del programa 

import csv
import matplotlib.pyplot as plt 
import numpy as np
import os
import pandas as pd

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
#los limites los saque de la tesis de edu, el resto elegi yo

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
                    aux = 1+3*w*c**2
                    A = (aux**(aux/(2*w*c**2)))/(4*w**(3/2)*c**3)
                    rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))
                    P = w*(c**2)*rho
        
                    dot_M_b = m_0*((2*a_b*t)/(3*(1-w)*T_b))*(1+(t/T_b)**2)**((-2-3*w)/(3*(1-w))) 
                    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P) 
                    if M > 10**17:
                        A_M = (5.3*10**25)        # en unidades de g^3 s^-1
                        dot_M_rh = - A_M / M**2
                    else: 
                        A_M = (7.8*10**26)       # en unidades de g^3 s^-1
                        dot_M_rh = - A_M / M**2
                else:
                    w = 1/3   
                    aux = 1+3*w*c**2
                    A = (aux**(aux/(2*w*c**2)))/(4*w**(3/2)*c**3)
                    rho = rho_0*a_b*(1+(t/T_b)**2)**(1/(3-3*w))  
                    P = w*(c**2)*rho                             
        
                    dot_M_b = m_0*((2*a_b*t)/(3*(1-w)*T_b))*(1+(t/T_b)**2)**((-2-3*w)/(3*(1-w))) 
                    dot_M_ac = 4*pi*A*G**2*c**(-5)*M**2*(c**2*rho-P) 
                    if M > 10**17:
                        A_M = (5.3*10**25)       # en unidades de g^3 s^-1
                        dot_M_rh = - A_M / M**2
                    else: 
                        A_M = (7.8*10**26)       # en unidades de g^3 s^-1
                        dot_M_rh = - A_M / M**2
        
            
                datos = guardar_datos_variacion(t, dot_M_b, dot_M_ac, dot_M_rh, nombre)
                
#%%
           
ruta = "/home/iara/faq/tesis/spyder"
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
 
    plt.plot(t, dot_M_b, label = "Dinámica espacio-tiempo con m_0 = {:.1e}".format(m_0), color = 'turquoise')
    plt.plot(t, dot_M_ac, label = "Acrecion", color='darkorange')
    plt.plot(t, dot_M_rh , label = "Radiacion", color = 'deeppink')
    plt.axvline(x = T_b, color = 'mediumpurple', label = 'escala del bounce')
    plt.legend()
    plt.title('Variacion considerando una masa M= {:.1e}'.format(M))
    plt.show()

    
#%%

