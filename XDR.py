#PROCESAMIENTO DE DATOS CRUDOS DE XDR - IDENTIFICACION Y CUANTIFICACION DE SEÑALES DE PHA
#Script creado por ANDREA TREJO y LUIS MARIO HERNANDEZ-SOTO 05/2025

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.signal import peak_widths
import math
from scipy.signal import find_peaks

# Cargar archivo
XDR = "8NaCl6AC_no_tab.tsv"  # Reemplaza con tu archivo
data = pd.read_csv(XDR, sep='\t')
#print(data.head())

# Extraer columnas
theta = data.iloc[:, 0].values  # Ángulo 2θ
intensidad = data.iloc[:, 1].values  # Intensive

# Graficar difractograma
plt.figure(figsize=(10, 6))
plt.plot(theta, intensidad, label='Difractograma')
plt.xlabel('2θ (°)')
plt.ylabel('Intensidad (u.a.)')
plt.title('Difracción de rayos X - PHA')
plt.legend()
plt.grid(True)
plt.savefig("difractograma.png", dpi=600)
plt.show()

#Paso 2: Encontrar los picos principales
# Encontrar picos
#picos, _ = find_peaks(intensidad, height=100, distance=10)  # Ajusta 'height' y 'distance' según tu muestra
# Umbral dinámico: percentil 75
umbral_pico = np.percentile(intensidad, 75)
picos, _ = find_peaks(intensidad, height=umbral_pico, distance=10)

# Mostrar picos detectados
for i in picos:
    print(f"Pico en 2θ = {theta[i]:.2f}° con intensidad = {intensidad[i]:.2f}")

# Paso 3: Calcular el índice de cristalinidad Método de dos áreas: comparar área de
# picos vs.fondo amorfo

# Área total bajo la curva
area_total = simpson(intensidad, theta)

# Umbral para separar cristalino y amorfo
umbral = np.percentile(intensidad, 50)

# Intensidad cristalina
intensidad_cristalina = np.where(intensidad > umbral, intensidad, 0)

# Área cristalina
area_cristalina = simpson(intensidad_cristalina, theta)

# Índice de cristalinidad
cristalinidad = (area_cristalina / area_total) * 100
print(f"Índice de cristalinidad: {cristalinidad:.2f}%")

# Paso 4: Calcular el tamaño de cristalita(ecuación de Scherrer)

# Parámetros para la ecuación de Scherrer
K = 0.9
lambda_radiacion = 1.5406  # Cu Kα en Å

# Tomar el primer pico fuerte para ejemplo
pico_principal = picos[0]
FWHM = peak_widths(intensidad, [pico_principal], rel_height=0.5)[0][0]
beta_rad = (FWHM * (np.pi / 180))  # Convertir a radianes
theta_rad = (theta[pico_principal] / 2) * (np.pi / 180)

# Cálculo del tamaño de cristalita
L = (K * lambda_radiacion) / (beta_rad * np.cos(theta_rad))
print(f"Tamaño de cristalita (Scherrer): {L:.2f} nm")

# Paso 5: Comparar con picos de PHB / PHBV
# Picos conocidos de PHB
picos_PHB = [13.5, 16.9, 22.5, 25.6]

# Comparar
print("\nComparación con picos PHB:")
for pico_ref in picos_PHB:
    diferencias = np.abs(theta[picos] - pico_ref)
    idx_min = np.argmin(diferencias)
    pico_encontrado = theta[picos[idx_min]]
    print(f"Pico teórico: {pico_ref:.1f}° — Pico encontrado: {pico_encontrado:.2f}°")

with open("resultados.txt", "w") as f:
    f.write("Resultados del análisis de difracción:\n\n")

    for i in picos:
        f.write(f"Pico en 2θ = {theta[i]:.2f}° con intensidad = {intensidad[i]:.2f}\n")

    f.write(f"\nÍndice de cristalinidad: {cristalinidad:.2f}%\n")
    f.write(f"Tamaño de cristalita (Scherrer): {L:.2f} nm\n")

    f.write("\nComparación con picos PHB:\n")
    for pico_ref in picos_PHB:
        diferencias = np.abs(theta[picos] - pico_ref)
        idx_min = np.argmin(diferencias)
        pico_encontrado = theta[picos[idx_min]]
        f.write(f"Pico teórico: {pico_ref:.1f}° — Pico encontrado: {pico_encontrado:.2f}°\n")
