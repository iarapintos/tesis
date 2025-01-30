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
eta_0 = 2.96e12 #son aproximadamente los 47k anios donde domina la rad
eta_b = 21969.3   #escala del bounce
a_b = 7.41155e-9 #factor de normalizacion del factor de escala
xini = eta_0/eta_b  #cte inicial en funcion de x

#definimos cte de normalizacion 
rho_rad_ini = (3*c**2/(8*pi*G))*(1/(a_b*eta_b)**2)*xini**2*(1+xini**2)**(-3)  #dens de rad inicial [erg/cm^3]
#beta = 1e-15
#beta_obs = 0.99999999999999 #el dato observacional 
#beta = beta_obs/(1-beta_obs)
beta =1e5
rho_bh_ini = beta*rho_rad_ini   #con el c^2 para que de las energias 

alfa = 2.5
M_min = m_sol*1e4
M_max = m_sol*1e6

cteA = (rho_bh_ini*(2-alfa))/(c**2*(M_max**(2-alfa)-M_min**(2-alfa)))

#intervalo de tiempo ,log
#ta_int= np.logspace(8,13,1000)
eta_pos = np.logspace(0, 13, 10000)
eta_neg = -eta_pos[::-1] 
eta_int = np.concatenate((eta_neg,eta_pos))

#aca viene la discretizacion

rho_bh = np.zeros_like(eta_int) 

Nm = 100
Masas = np.logspace(np.log10(M_min), np.log10(M_max),Nm)

hlog = (M_max/M_min)**(1/(Nm-1))

i = 0
for eta in eta_int:
    a = a_b*(1+(eta/eta_b)**2)**(1/2) #factor de escala
    for M in Masas:
        M_i = M*a
        N_i = cteA*M**(-alfa)*a**(-3)   #densidad inicial de bhs en el bloque i

        delta = M*(hlog-1)

        rho_bh[i] = rho_bh[i]+ N_i*M_i*c**2*delta
        
    i = i+1

# plot results
plt.plot(eta_int, rho_bh, color = 'deeppink')
plt.xlabel(r'$\eta$')
#plt.xlim([1e8, 1e15])
plt.xscale('symlog')
plt.ylabel(r'$\rho_{bh}$')
plt.yscale('log')
plt.axhline(y=rho_bh_ini, color='turquoise', linestyle='-')
plt.show()

#ahora vamos con la densidad del fluido cosmologico 

ctes = (3*c**2)/(8*pi*G*a_b**2*eta_b**2)

rho_cf = ctes*(eta_int/eta_b)**2*(1+(eta_int/eta_b)**2)**(-3)-rho_bh
rho_cf_solo = ctes*(eta_int/eta_b)**2*(1+(eta_int/eta_b)**2)**(-3)

# plot results
plt.plot(eta_int, rho_cf, color = 'deeppink', label=r'$\rho_{cf}$')
plt.plot(eta_int, rho_cf_solo, color = 'darkorange', label=r'$\rho_{cf}$ sin $\rho_{bh}$')
plt.plot(eta_int, rho_bh, color = 'turquoise', label=r'$\rho_{bh}$')
plt.xlabel(r'$\eta$')
#plt.xlim([1e8, 1e15])
plt.xscale('symlog')
plt.ylabel(r'$\rho_{cf}$')
plt.yscale('log')
plt.legend()
plt.show()


#omega cf 

#primero defino la derivda de la densidad del fluido cosmologico

def der_rho(rho, etavar, k):
    drho = (rho[k+1] - rho[k])/(etavar[k+1]-etavar[k])
    return drho

der_rho_cf = np.zeros_like(eta_int) 
for k in range(len(eta_int)-1):
    der_rho_cf[k]= der_rho(rho_cf, eta_int, k)


#defino Q

Q_inter = np.zeros_like(eta_int) 

a_vector = a_b*(1+(eta_int/eta_b)**2)**(1/2) #factor de escala definido como vector

j = 0
for eta in eta_int:
    a = a_b*(1+(eta/eta_b)**2)**(1/2) #factor de escala
    for M in Masas:
        M_i = M*a
        N_0 = cteA*M**(-alfa)

        dNdeta = N_0/(eta_int[j+1]-eta_int[j])*(a_vector[j+1]**(-3)-a_vector[j]**(-3))

        delta = M*(hlog-1)

        Q_inter[j] = Q_inter[j]+ dNdeta*M_i*c**2*delta
        
    j = j+1 if j != 2000 else 1999

omega_cf = (1/3)*(eta_b**2/eta_int)*(1+(eta_int/eta_b)**2)*(1/rho_cf)*(-Q_inter-der_rho_cf)-1

# plot results
plt.plot(eta_int, omega_cf, color = 'deeppink')
plt.xlabel(r'$\eta$')
#plt.xlim([1e8, 1e15])
plt.xscale('symlog')
plt.ylabel(r'$\omega_{cf}$')
plt.yscale('symlog', linthresh=0.000000001)
#plt.axhline(y=0, color='turquoise', linestyle='-')
plt.show()
