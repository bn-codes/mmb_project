# Načtení knihoven
from math import factorial, exp
from decimal import Decimal, getcontext
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

# Nastavení požadované přesnosti čísel
getcontext().prec = 50

# Definice zadaných parametrů

# Intervaly věku v letech
age_groups = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85,
              90, 95]

# Průměrný počet mutací pro každou buňku
generation = [45, 46, 47, 47.5, 48, 48.5, 49, 49.5, 50, 50.5, 51, 51.5, 52, 
              52.5, 53, 53.5, 54, 54.5, 55]

# Počet obnov buněk během jednoho roku
parameter_n = [4.2e+13, 4.2e+13, 4.2e+13, 4.2e+13, 4.2e+13, 4.2e+13, 4.2e+13,
               3.15e+13, 2.36e+13, 1.77e+13, 1.33e+13, 9.97e+12, 7.48e+12, 5.61e+12,
               4.2e+12, 2.73e+12, 1.78e+12, 1.07e+12, 5.33e+11]
parameter_n_decimal = [Decimal(n) for n in parameter_n]

# Pravděpodobnost, že se zdravá buňka změní v rakovinnou v průběhu jednoho dělení
p_constant = Decimal(2.38e-18)
p_constant_list = [p_constant] * len(parameter_n)

# Pozorované pravděpodobnosti výskytu rakoviny za 5 let podle Cancerstats UK
cancerstats = [0.1028, 0.0553, 0.0633, 0.1023, 0.1643, 0.3003, 0.4533, 0.6380, 0.9550, 1.5588, 
               2.3953, 3.5565, 5.3138, 7.5760, 9.5098, 11.8208, 13.0510, 14.2038, 13.3100]

# Práh počtu akumulovaných mutací, při jehož překročení se buňka stává rakovinnou
q = 118

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
   
    Parametry:
        n (Decimal): Počet obnov buněk.
        pc (Decimal): Pravděpodobnost vzniku rakoviny během jednoho dělení.
        pa (float): Pravděpodobnost vzniku rakoviny kvůli akumulaci škodlivých mutací.
   
   Návratová hodnota:
       float: Pravděpodobnost zachování zdraví.
   """
   
    result = 1 / (exp(float(n * (pc + pa))))
    return result


# Dopočítání ostatních proměnných potřebných pro model

# Pravděpodobnost vzniku rakoviny kvůli akumulaci škodlivých mutací
p_accumulative = [compute_pa(q, lambd) for lambd in generation]

# Součet pravděpodobnosti vzniku rakoviny během jednoho dělení a kvůli akumulaci škodlivých mutací
p_all = [x + p_constant for x in p_accumulative]


# Výpočet p_cancer_5_years pro zadané nekonstantní parametry n

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


# Výpočet p_cancer_5_years pro konstantní parameter_n = 4.2e+13
parameter_n_constant = [Decimal(4.2e+13)] * len(parameter_n)

# Konečná pravděpodobnost vzniku rakoviny
lambda_np_constant = [n * p for n, p in zip(parameter_n_constant, p_all)]

# Exponent konečné pravděpodobnosti vzniku rakoviny
e_lambda_constant = [exp(float(x)) for x in lambda_np_constant]

# Pravděpodobnost, že všechny buňky jsou zdravé
p_0_constant = [compute_phealth(n, pc, pa) for n, pc, pa in zip(parameter_n_constant, 
                                                                p_constant_list, p_accumulative)]

# Pravděpodobnost vzniku rakoviny v průběhu jednoho roku života
p_cancer_year_constant = [1 - x for x in p_0_constant]

# Pravděpodobnost (v %) vzniku rakoviny v průběhu pěti let života
p_cancer_5_years_constant = [100 * 5 * x for x in p_cancer_year_constant]

# Funkce pro výpočet R-squared
def calculate_r2(true_values, predicted_values):
    
    """
    Funkce pro výpočet koeficientu determinace (R^2).

    Parametry:
        true_values (array-like): Skutečné hodnoty.
        predicted_values (array-like): Predikované hodnoty.

    Návratová hodnota:
        float: Koeficient determinace (R^2), který měří, jak dobře predikované hodnoty odpovídají skutečným hodnotám.
    """
    
    return r2_score(true_values, predicted_values)

# Výpočet R-squared pro dva různé modely
r2_model_with_turnover = calculate_r2(cancerstats, p_cancer_5_years)
r2_model_without_turnover = calculate_r2(cancerstats, p_cancer_5_years_constant)

# Vytvoření grafu
plt.figure(figsize=(10, 6))

# Vykreslení p_cancer_5_years pro zadané nekonstantní parametry n
plt.plot(age_groups, p_cancer_5_years, label=f'Model with turnover reduction  R²={r2_model_with_turnover:.3f}', color='red', marker='o')

# Vykreslení p_cancer_5_years pro konstantní parametr n
plt.plot(age_groups, p_cancer_5_years_constant, label=f'Model without turnover reduction  R²={r2_model_without_turnover:.3f}', color='blue', marker='o')

# Vykreslení cancerstats
plt.plot(age_groups, cancerstats, label='Cancerstats', linestyle='--', color='black', marker='^')

# Nastavení osy y
plt.yscale('log')
plt.ylim(0.02, 110)
plt.yticks([0.1, 1, 10, 100], ['0.1', '1', '10', '100'])

# Přidání názvů a legendy
plt.xlabel('Age')
plt.ylabel('Cancer rates/5 years')
plt.title('Cancer Probability with/out turnover reduction')
plt.legend()

# Zobrazení grafu
plt.show()
