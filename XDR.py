#PROCESAMIENTO DE DATOS CRUDOS DE XDR - IDENTIFICACION Y CUANTIFICACION DE SEÑALES DE PHA
#Script creado por ANDREA TREJO y LUIS MARIO HERNANDEZ-SOTO 05/2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import os

# Constantes
K_SCHERRER = 0.9
LONGITUD_ONDA = 1.5406  # Å (Cu Kα)
PHB_PEAKS = [13.5, 16.9, 22.5, 25.6]

# Archivos
archivos = [f for f in os.listdir() if f.endswith('.tsv')]
salida_txt = open("resultados_cristalinidad.txt", "w")

for archivo in archivos:
    df = pd.read_csv(archivo, sep="\t")
    angulo = df.iloc[:, 0].values
    intensidad = df.iloc[:, 1].values

    # Detección de picos solo con distancia mínima
    picos, _ = find_peaks(intensidad, distance=10)

    # Calcular índice de cristalinidad
    if len(picos) > 0:
        intensidad_total = np.sum(intensidad)
        intensidad_cristalina = np.sum(intensidad[picos])
        indice_cristalinidad = (intensidad_cristalina / intensidad_total) * 100
    else:
        indice_cristalinidad = np.nan

    # Calcular tamaño de cristalita
    if len(picos) > 0:
        pico_max = picos[np.argmax(intensidad[picos])]
        fwhm = 0.2  # Valor fijo o estimado real si tienes los datos
        theta_rad = np.deg2rad(angulo[pico_max] / 2)
        tamano_cristalita = (K_SCHERRER * LONGITUD_ONDA) / (fwhm * np.cos(theta_rad))
    else:
        tamano_cristalita = np.nan

    # Comparación con picos de PHB
    coincidencias = []
    for pico_teorico in PHB_PEAKS:
        if len(picos) > 0:
            diferencias = np.abs(angulo[picos] - pico_teorico)
            if np.min(diferencias) < 0.3:
                coincidencias.append(f"Pico teórico: {pico_teorico}° — Pico encontrado: {angulo[picos][np.argmin(diferencias)]:.2f}°")
            else:
                coincidencias.append(f"Pico teórico: {pico_teorico}° — No se detectaron picos")
        else:
            coincidencias.append(f"Pico teórico: {pico_teorico}° — No se detectaron picos")

    # Imprimir resultados al archivo
    salida_txt.write(f"Archivo: {archivo}\n")
    if len(picos) > 0:
        for i in picos:
            salida_txt.write(f"  Pico en 2θ = {angulo[i]:.2f}° con intensidad = {intensidad[i]:.2f}\n")
    salida_txt.write(f"  Índice de cristalinidad: {indice_cristalinidad:.2f}%\n" if not np.isnan(indice_cristalinidad) else "  Índice de cristalinidad: nan%\n")
    salida_txt.write(f"  Tamaño de cristalita (Scherrer): {tamano_cristalita:.2f} nm\n" if not np.isnan(tamano_cristalita) else "  No se detectó pico principal para calcular Scherrer.\n")
    salida_txt.write("  Comparación con picos PHB:\n")
    for c in coincidencias:
        salida_txt.write(f"    {c}\n")
    salida_txt.write("\n" + "-"*60 + "\n\n")

    # Graficar
    plt.figure(figsize=(10, 5))
    plt.plot(angulo, intensidad, label="Difractograma")
    if len(picos) > 0:
        plt.plot(angulo[picos], intensidad[picos], "rx", label="Picos detectados")
    plt.xlabel("2θ (°)")
    plt.ylabel("Intensidad (u.a.)")
    plt.title(f"Difractograma: {archivo}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{archivo}_difractograma.png")
    plt.close()

salida_txt.close()
print("Análisis completo. Resultados en 'resultados_cristalinidad.txt'")
