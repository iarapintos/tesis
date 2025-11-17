import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction

# Leer el CSV especificando que los campos están entre comillas
data = pd.read_csv("omega_sin_bh.csv", header=None, names=["x", "omega"], quotechar='"')

# Convertir las fracciones a floats con cuidado (esto puede tardar si no se limita)
def frac_to_float(s):
    try:
        return float(Fraction(s))
    except:
        return np.nan

# Aplica la conversión de forma más controlada
data["x"] = data["x"].apply(frac_to_float)
data["omega"] = data["omega"].apply(frac_to_float)

# Eliminar filas con errores de conversión
data = data.dropna()

# Si querés ver solo una parte para evitar cuelgues
# data = data.iloc[:5000]  # descomenta para probar con menos puntos

# Graficar
plt.plot(data["x"], data["omega"], color='mediumseagreen')
plt.xlabel(r'$\eta/\eta_b$')
plt.ylabel(r'$\omega_{FC}$')
plt.axhline(1/3, color='orchid', linestyle='--', linewidth=1)
plt.xlim([-5,5])
plt.ylim([-1e1,2])

plt.show()