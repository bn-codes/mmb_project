# Načtení knihoven
import numpy as np
import matplotlib.pyplot as plt

# Definice konstanty p (pravděpodobnost vzniku rakoviny)
p = 1e-15

# Definice funkce P-cancer jako pravděpodobnosti vzniku rakoviny 
# v závislosti na počtu dělení buněk
def P_cancer(x, p):
    
    """
    Funkce pro výpočet pravděpodobnosti vzniku rakoviny.
    
    Parametry:
        x (array-like): Počet dělení buněk.
        p (float): Konstanta pravděpodobnosti.
    
    Návratová hodnota:
        array-like: Pravděpodobnost vzniku rakoviny.
    """
    
    return 1 - (1 / np.exp(2**x * p))

# Určení hodnot na ose x
x_values = np.linspace(0, 80, 1000)

# Výpočet pravděpodobnosti rakoviny pro dané hodnoty x
y_values = P_cancer(x_values, p)

# Vytvoření grafu
plt.figure(figsize=(10, 6))

# Vykreslení výsledků
plt.plot(x_values, y_values, label=rf'$P_{{cancer}} = 1 - \frac{{1}}{{e^{{2^x \cdot p}}}}, \ p = {p}$', color='r')

# Přidání názvů a legendy
plt.title('P(rakovinné) jako funkce počtu dělení buněk' )
plt.xlabel('Počet dělení buněk')
plt.ylabel('P_{rakovinné}')
plt.legend()

# Zobrazení grafu
plt.grid(True)
plt.show()
