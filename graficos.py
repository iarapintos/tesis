import pandas as pd
import matplotlib.pyplot as plt

from fractions import Fraction


# DENSIDAD DE AGUEJOS NEGROS
data = pd.read_csv("rhobh.csv", header=None, names=["x", "rhobh"])

plt.plot(data["x"], data["rhobh"], color = 'deeppink')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\rho_{AN}/\rho_{AN}(\eta_{ci}/\eta_b)$')
plt.yscale('log')
plt.xlim([-5,5])
plt.show()

#DENSIDAD DE ENERGIA DEL FLUDO COSMOLOGICO
data = pd.read_csv("rhocf.csv", header=None, names=["x", "rhocf"])

plt.plot(data["x"], data["rhocf"], color = 'turquoise')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\rho_{FC}/\rho_{FC}(\eta_{ci}/\eta_b)$')
plt.yscale('log')
plt.show()


#DENSIDAD CF CERCA 



#DENSIDAD CF LINDA 

data = pd.read_csv("rhocfcerca.csv", header=None, names=["x", "rhocfcerca"])
data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

data_sinbh = pd.read_csv("rhocf_cerca_sin_gasBH.csv", header=None, names=["x", "rhocf_sinBH"])
data_sinbh["x"] = data_sinbh["x"].apply(lambda x: float(Fraction(x)))

fig, ax = plt.subplots()

plt.plot(data["x"], data["rhocfcerca"], color = 'turquoise', label=r'$\rho_{cf}$')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\rho_{cf}/\rho_{cf}(\eta_{ci}/\eta_b)$',  labelpad=-3)
plt.yscale('symlog')
plt.xscale('symlog', linthresh=0.00000000001)
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.xlim([-5e-9, 5e-9])

intervalo = [ 1e2, 1e6, 1e10, 1e14, 1e18]
int_neg = -np.array(intervalo)[::-1]

ticks_y = np.concatenate((int_neg, [0], intervalo))
ax.set_yticks(ticks_y)

plt.show()


#OMEGA 


data = pd.read_csv("omega.csv", header=None, names=["x", "omega"])

plt.plot(data["x"], data["omega"], color = 'mediumseagreen')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\omega_{FC}$')
plt.axhline(1/3, color='orchid', linestyle='--', linewidth=1)
plt.xlim([-5, 5])
plt.ylim([-2, 1])
plt.show()

#OMEGA CERCA Y DE A DOS 

data = pd.read_csv("omeganeg.csv", header=None, names=["x", "omega"])

data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

plt.plot(data["x"], data["omega"], color = 'mediumseagreen')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\omega_{cf}$')
plt.xscale('symlog', linthresh=0.00000000000001)
plt.yscale('symlog', linthresh=0.00000000000001)
plt.show()



data = pd.read_csv("omegapos.csv", header=None, names=["x", "omega"])

data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

plt.plot(data["x"], data["omega"], color = 'mediumseagreen')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\omega_{cf}$')
plt.xscale('symlog', linthresh=0.000000000000000001)
plt.yscale('log')
plt.show()

#CONDICION NEC

fig, ax = plt.subplots()

data = pd.read_csv("NEC.csv", header=None, names=["x", "NEC"])
data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

plt.plot(data["x"], data["NEC"], color = 'darkorange')
plt.xlabel(r'$\eta/\eta_b$')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.ylabel(r'$\rho_{FC}+P_{FC}~/c^2$', labelpad=-3)
plt.xscale('symlog', linthresh=0.0001)
plt.yscale('symlog')
plt.ylim([-1e38,1e35])

intervalo = [ 1e5, 1e10, 1e15, 1e20, 1e25, 1e30, 1e35]
int_neg =  [-1e5, -1e10, -1e15, -1e20, -1e25, -1e30, -1e35, -1e40]

ticks_y = np.concatenate((int_neg, [0], intervalo))
ax.set_yticks(ticks_y)

plt.show()


#CONDICION SEC

data = pd.read_csv("SEC.csv", header=None, names=["x", "SEC"])
data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

plt.plot(data["x"], data["SEC"], color = 'mediumorchid')
plt.xlabel(r'$\eta/\eta_b$')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.ylabel('Condición SEC')
plt.xscale('symlog', linthresh=0.001)
plt.yscale('symlog')
plt.ylim([-1e40,1e35])
plt.show()


from fractions import Fraction

data = pd.read_csv("omega_sin_gasBH.csv", header=None, names=["x", "SEC"])
data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

plt.plot(data["x"], data["SEC"], color = 'mediumseagreen')
plt.xlabel(r'$\eta/\eta_b$')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.ylabel(r'$\omega_{CF}$')
plt.xlim([-3,3])
plt.ylim([-5, 2])
plt.show()


# AMBOS DENSIDAD CF MONOENERGETICA LEJOS
data = pd.read_csv("rhocf.csv", header=None, names=["x", "rhocf"])

data_sinbh = pd.read_csv("rhocf_sin_gasBH.csv", header=None, names=["x", "rhocf_sinBH"])

plt.plot(data["x"], data["rhocf"], color = 'turquoise', label=r'$\rho_{cf}$')
plt.plot(data_sinbh['x'], data_sinbh["rhocf_sinBH"], color= 'darkorange', label=r'$\rho_{cf}$ en ausencia de BH')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\rho_{cf}$ [erg/cm$^3$]')
plt.yscale('log')
plt.legend()
plt.show()


#AMBAS DENSIDADES DE ENERGIA CF CERCA DEL BOUNCE 


data = pd.read_csv("rhocfcerca.csv", header=None, names=["x", "rhocfcerca"])
data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

data_sinbh = pd.read_csv("rhocf_cerca_sin_gasBH.csv", header=None, names=["x", "rhocf_sinBH"])
data_sinbh["x"] = data_sinbh["x"].apply(lambda x: float(Fraction(x)))

fig, ax = plt.subplots()

plt.plot(data["x"], data["rhocfcerca"], color = 'turquoise', label=r'$\rho_{cf}$')
plt.plot(data_sinbh['x'], data_sinbh["rhocf_sinBH"], color= 'darkorange', label=r'$\rho_{cf}$ en ausencia de BH')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\rho_{cf}/\rho_{cf}(\eta_{ci}/\eta_b)$',  labelpad=-3)
plt.yscale('symlog', linthresh=0.00000000001)
plt.xscale('symlog', linthresh=0.00000000001)
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.xlim([-5e-9, 5e-9])
plt.legend(loc='lower center', bbox_to_anchor=(0.5, 0.15))

intervalo = [ 1e-7, 1e-3, 1e1, 1e5, 1e9, 1e13, 1e17, 1e21]
int_neg = -np.array(intervalo)[::-1]

ticks_y = np.concatenate((int_neg, [0], intervalo))
ax.set_yticks(ticks_y)

plt.show()


# OMEGAS CERCA DEL BOUNCE 


fig, ax = plt.subplots()

data = pd.read_csv("omeganeg.csv", header=None, names=["x", "omega"])

data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

plt.plot(data["x"], data["omega"], color = 'mediumseagreen')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\omega_{FC}$')
plt.xscale('symlog', linthresh=0.00000000000001)
plt.yscale('symlog')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.xlim(-1e-6, 1e-6)

intervalo = [ 1e3, 1e6, 1e9, 1e12, 1e15, 1e18, 1e21]
int_neg = -np.array(intervalo)[::-1]

ticks_y = np.concatenate((int_neg, [0], intervalo))
ax.set_yticks(ticks_y)

plt.xticks([-1e-7, -1e-9, -1e-11, -1e-13, 0, 1e-13, 1e-11, 1e-9, 1e-7], [r'$-10^{-7}$', r'$-10^{-9}$', r'$-10^{-11}$', r'$-10^{-13}$', '0',  r'$10^{-13}$', r'$10^{-11}$', r'$10^{-9}$', r'$10^{-7}$']
)

plt.show()



#CONDICION SEC PIOLA

fig, ax = plt.subplots()

data = pd.read_csv("SEC.csv", header=None, names=["x", "SEC"])
data["x"] = data["x"].apply(lambda x: float(Fraction(x)))

plt.plot(data["x"], data["SEC"], color = 'mediumorchid')
plt.xlabel(r'$\eta/\eta_b$')
plt.axhline(0, color='gray', linestyle='--', linewidth=1)
plt.ylabel(r'$\rho_{FC}+3P_{FC}~/c^2$', labelpad=-3)
plt.xscale('symlog', linthresh=0.0001)
plt.yscale('symlog')
plt.ylim([-1e38,1e35])

intervalo = [ 1e5, 1e10, 1e15, 1e20, 1e25, 1e30, 1e35]
int_neg =  [-1e5, -1e10, -1e15, -1e20, -1e25, -1e30, -1e35, -1e40]

ticks_y = np.concatenate((int_neg, [0], intervalo))
ax.set_yticks(ticks_y)

plt.show()

#factor de escala 


eta_b = 21969.3   #escala del bounce
a_b = 7.41155e-9 #factor de normalizacion del factor de escala

eta = np.linspace(-8e5, 8e5, 2000)

a = a_b*(1+(eta/eta_b)**2)**(1/2) #factor de escala

fig, ax = plt.subplots()

plt.plot(eta, a, color = 'turquoise')
plt.xlabel(r'$\eta$ [s]')
plt.ylabel(r'$a$')
plt.yscale('log')
#plt.xscale('symlog')


plt.xticks([-8e5, -6e5, -4e5, -2e5, 0, 2e5, 4e5, 6e5, 8e5], [r'$-8 \times 10^5$', r'$-6 \times 10^5$', r'$-4 \times 10^5$',  r'$-2 \times 10^5$', '0', r'$2 \times 10^5$', r'$4 \times 10^5$',  r'$6 \times 10^5$', r'$8 \times 10^5$'],fontsize=8  # Ajusta este valor (prueba con 8, 6, etc.) 
           )

plt.show()

