# Načtení knihoven
from math import factorial, exp
from decimal import Decimal, getcontext
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

# Nastavení požadované přesnosti čísel
getcontext().prec = 50

# Definice zadaných parametrů

# Intervaly věku v letech
age_groups = [5, 10, 15, 20, 25, 30, 35]

# Průměrný počet mutací pro každou buňku
generation = [45, 46, 47, 47.5, 48, 48.5, 49]

# Počet obnov buněk během jednoho roku
parameter_n = [4.2e+13, 4.2e+13, 4.2e+13, 4.2e+13, 4.2e+13, 4.2e+13, 4.2e+13]
parameter_n_decimal = [Decimal(n) for n in parameter_n]

# Pravděpodobnost, že se zdravá buňka změní v rakovinnou v průběhu jednoho dělení
p_constant = Decimal(2.38e-18)
p_constant_list = [p_constant] * len(parameter_n)

# Pozorované pravděpodobnosti výskytu rakoviny za 5 let podle Cancerstats UK
cancerstats = [0.1028, 0.0553, 0.0633, 0.1023, 0.1643, 0.3003, 0.4533]

# Různé prahy počtu akumulovaných mutací, při jehož překročení se buňka stává rakovinnou
q_list = [117, 118, 119]

# Definice funkce pro výpočet p_accumulative
# Pravděpodobnost, že se zdravá buňka změní v rakovinnou po akumulaci mutací ve vícero děleních
def compute_pa(q, lambd):
    
    """
    Funkce pro výpočet pravděpodobnosti vzniku rakoviny při akumulaci škodlivých mutací.
    
    Parametry:
        q (int): Práh mutací.
        lambd (Decimal): Parametr rozdělení.
    
    Návratová hodnota:
        float: Pravděpodobnost vzniku rakoviny kvůli akumulaci škodlivých mutací.
    """
    
    lambd = Decimal(lambd) 
    summation = sum((lambd**Decimal(k)) / Decimal(factorial(k)) for k in range(q + 1))
    pa = 1 - summation / lambd.exp()
    return pa


# Definice funkce pro výpočet p_health (pravděpodobnost, že všechny buňky jsou zdravé)
def compute_phealth(n, pc, pa):
    
    """
    Funkce pro výpočet pravděpodobnosti zachování zdraví.
   # Součet pravděpodobnosti vzniku rakoviny během jednoho dělení a kvůli akumulaci škodlivých mutací
    Parametry:
        n (Decimal): Počet obnov buněk.
        pc (Decimal): Pravděpodobnost vzniku rakoviny během jednoho dělení.
        pa (float): Pravděpodobnost vzniku rakoviny kvůli akumulaci škodlivých mutací.
   
   Návratová hodnota:
       float: Pravděpodobnost zachování zdraví.
   """
   
    result = 1 / (exp(float(n * (pc + pa))))
    return result

# Prázdný slovník pro uložení výsledků pro různé hodnoty q
results = {}

# Pro všechny možné hodnoty prahů
for q in q_list:
    
    # Dopočítání ostatních proměnných potřebných pro model
    
    # Pravděpodobnost vzniku rakoviny kvůli akumulaci škodlivých mutací
    p_accumulative = [compute_pa(q, lambd) for lambd in generation]
    
    # Součet pravděpodobnosti vzniku rakoviny během jednoho dělení a kvůli akumulaci škodlivých mutací
    p_all = [x + p_constant for x in p_accumulative]
    
    # Konečná pravděpodobnost vzniku rakoviny
    lambda_np = [n * p for n, p in zip(parameter_n_decimal, p_all)]
    
    # Exponent konečné pravděpodobnosti vzniku rakoviny
    e_lambda = [exp(float(x)) for x in lambda_np]
    
    # Pravděpodobnost, že všechny buňky jsou zdravé
    p_0 = [compute_phealth(n, pc, pa) for n, pc, pa in zip(parameter_n_decimal, 
                                                           p_constant_list, p_accumulative)]
    
    # Pravděpodobnost vzniku rakoviny v průběhu jednoho roku života
    p_cancer_year = [1 - x for x in p_0]
    
    # Pravděpodobnost (v %) vzniku rakoviny v průběhu pěti let života
    p_cancer_5_years = [100 * 5 * x for x in p_cancer_year]
    
    # Uložení výsledků pro daný práh
    results[q] = p_cancer_5_years
    
    
# Vytvoření grafu
plt.figure(figsize=(10, 6))

# Definice markerů a barev
markers = ['s', 'o', '^', 'v']
colors = ['blue', 'red', 'lightblue', 'black']

# Vykreslení výsledků pro každé q, včetně výpočtu R-squared
for i, q in enumerate(q_list):
    
    # Výpočet R-squared mezi modelovými výsledky a cancerstats
    r_squared = r2_score(cancerstats, results[q])
    
    # Přidání R-squared do labels
    plt.plot(age_groups, results[q], 
             label=f'Model q = {q}, R² = {r_squared:.3f}', 
             marker=markers[i % len(markers)], 
             color=colors[i % len(colors)])

# Přidání cancerstats do grafu
plt.plot(age_groups, cancerstats, label='Cancerstats', linestyle='--', color=colors[-1], marker=markers[-1])

# Přidání názvů a legendy
plt.xlabel('Age')
plt.ylabel('Cancer rates/5 years')
plt.title('Comparison of Cancer Probability for Different q Values')
plt.legend()

# Nastavení logaritmické osy y a omezení jejího rozsahu
plt.yscale('log')
plt.ylim(0.01, 1.2)

# Nastavení hodnot ticků pro logaritmickou osu y
plt.yticks([0.01, 0.1, 1], ['0.01', '0.1', '1'])

# Zobrazení grafu
plt.show()


