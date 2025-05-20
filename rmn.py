#PROCESAMIENTO DE DATOS CRUDOS RMN - IDENTIFICACION Y CUANTIFICACION DE SEÑALES DE PHA
#Script creado por ANDREA TREJO y LUIS MARIO HERNANDEZ-SOTO 05/2025

import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.integrate import trapezoid

#ruta a la carpeta donde estan los archivos
bruker_dir="."
dic, data=ng.bruker.read(bruker_dir)

#Apodización (line broadening) + Fourier Transform
lb= 1 #Hz
#print(dic.keys())
#print(dic['acqus'].keys())

fid = data * np.exp(-lb * np.arange(data.shape[0]) / dic['acqus']['SW_h'])

#Transformada de Fourier
spectrum = np.fft.fftshift(np.fft.fft(fid))
spectrum = np.abs(spectrum)

#Eje de ppm
sw = dic['acqus']['SW_h']       # ancho de espectro en Hz
o1 = dic['acqus']['O1']         # offset en Hz
bf1 = dic['acqus']['BF1']       # frecuencia espectrómetro en MHz
n_points = dic['acqus']['TD']   # puntos adquiridos
#hz_axis = np.linspace(-sw/2, sw/2, n_points)
#ppm_axis = (hz_axis + o1) / bf1
# Parámetros espectrales
sw = dic['acqus']['SW_h']       # ancho de espectro en Hz
o1 = dic['acqus']['O1']         # offset en Hz
bf1 = dic['acqus']['BF1']       # frecuencia espectrómetro en MHz
n_points = spectrum.size        # usa la longitud real del espectro

# Eje de frecuencias en Hz (centrado)
hz_axis = np.linspace(-sw/2, sw/2, n_points)
ppm_axis = (hz_axis + o1) / bf1  # conversión a ppm

# Graficar espectro
plt.figure(figsize=(10, 5))
plt.plot(ppm_axis, spectrum)
plt.xlim(7, 0)  # Escala típica 1H
plt.xlabel("Desplazamiento químico (ppm)")
plt.ylabel("Intensidad")
plt.title("Espectro 1H RMN de PHA")
plt.grid()
plt.show()

#Para cuantificar PHA
def integrar_region(ppm_axis, spectrum, rango_ppm):
    """Integra el área bajo la curva en un rango ppm dado."""
    mask = (ppm_axis >= min(rango_ppm)) & (ppm_axis <= max(rango_ppm))
    area = np.trapezoid(spectrum[mask], ppm_axis[mask])
    return area

def cuantificar_PHA(ppm_axis, spectrum, usar_CH=True):
    """
Cuantifica proporción molar 3HB y 3HV a partir del espectro ¹H RMN.
    
    usar_CH=True: integra las señales de CH (metino, δ≈5.2 ppm)
    usar_CH=False: integra CH₃ (δ≈1.2 ppm para 3HB, 0.9 ppm para 3HV)
    """
    if usar_CH:
        # Señales de protones CH (un protón por unidad monomérica)
        rango_3HB = (5.15, 5.30)  # ppm
        rango_3HV = (5.15, 5.30)  # se solapan, se usa CH₃ para distinguir mejor si hay mezcla
    else:
        # Señales de grupos metilo CH₃
        rango_3HB = (1.15, 1.35)  # ppm
        rango_3HV = (0.85, 1.05)  # ppm

    area_3HB = integrar_region(ppm_axis, spectrum, rango_3HB)
    area_3HV = integrar_region(ppm_axis, spectrum, rango_3HV)

    total = area_3HB + area_3HV
    porcentaje_3HB = (area_3HB / total) * 100 if total > 0 else 0
    porcentaje_3HV = (area_3HV / total) * 100 if total > 0 else 0

    print(f"Área 3HB: {area_3HB:.2f}")
    print(f"Área 3HV: {area_3HV:.2f}")
    print(f"→ Porcentaje 3HB: {porcentaje_3HB:.1f} %")
    print(f"→ Porcentaje 3HV: {porcentaje_3HV:.1f} %")

    return porcentaje_3HB, porcentaje_3HV

#Una vez que tenga el espectro procesado como: 
ppm_axis  # eje x (ppm)
spectrum  # intensidad (abs)
#Tengo que llamar:
cuantificar_PHA(ppm_axis, spectrum, usar_CH=False)

# Guardar el espectro en un archivo CSV
df = pd.DataFrame({'ppm': ppm_axis, 'intensidad': spectrum})
df.to_csv("espectro_PHA.csv", index=False)
print("✔ Espectro guardado como 'espectro_PHA.csv'")

# Guardar la figura como imagen PNG
plt.figure(figsize=(10, 5))
plt.plot(ppm_axis, spectrum)
plt.xlim(7, 0)
plt.xlabel("Desplazamiento químico (ppm)")
plt.ylabel("Intensidad")
plt.title("Espectro 1H RMN de PHA")
plt.grid()
plt.tight_layout()
plt.savefig("espectro_PHA.png", dpi=300)
print("✔ Imagen guardada como 'espectro_PHA.png'")

######
# === FFT y procesamiento ===
fid = data  # asegúrate que no esté recortado
n = fid.shape[0]
lb = 1
fid = fid * np.exp(-lb * np.arange(n) / dic['acqus']['SW_h'])
spectrum = np.fft.fftshift(np.fft.fft(fid))
spectrum = np.abs(spectrum)

# === Suavizado y corrección de línea base ===
from scipy.signal import savgol_filter
baseline = savgol_filter(spectrum, 101, 3)
spectrum_corrected = np.clip(spectrum - baseline, 0, None)
spectrum_smoothed = savgol_filter(spectrum_corrected, 31, 2)

# === Eje de desplazamiento químico ===
sw = dic['acqus']['SW_h']
o1 = dic['acqus']['O1']
bf1 = dic['acqus']['BF1']
ppm_axis = np.linspace(-sw/2, sw/2, n)  # usa el mismo n que usaste para el fid
ppm_axis = (ppm_axis + o1) / bf1

# === Verificación (opcional) ===
print("ppm_axis:", len(ppm_axis))
print("spectrum_smoothed:", len(spectrum_smoothed))



# SEÑALES EN ESPECTRO 1H RMN
print(f"Senales 1H")
# === Corrección de línea base (subtracción de mínimo local) ===
baseline = savgol_filter(spectrum, 101, 3)  # ventana 101, polinomio grado 3
#print(f"linea_base")

#print(f"linea_base")
# === Asignaciones esperadas para PHA (3HB y 3HV) ===
asignaciones = {
    'CH₃ (3HV)': 0.9,
    'CH₃ (3HB)': 1.2,
    'CH₂ (3HV)': 1.6,
    'CH₂ (3HB)': 2.5,
    'CH (3HB/3HV)': 5.2,
}

# === Gráfica del espectro con anotaciones ===
plt.figure(figsize=(12, 6))
plt.plot(ppm_axis, spectrum_smoothed, color='black', label='Espectro procesado')
plt.xlim(6, 0.5)
plt.xlabel("Desplazamiento químico (ppm)")
plt.ylabel("Intensidad (a.u.)")
plt.title("¹H RMN de PHA: Asignación de Señales")
plt.grid(True)

# Marcar señales esperadas
for nombre, ppm in asignaciones.items():
    plt.axvline(ppm, color='red', linestyle='--', alpha=0.6)
    plt.text(ppm + 0.05, max(spectrum_smoothed)*0.6, nombre,
             rotation=90, verticalalignment='center', color='red', fontsize=9)

plt.legend()
plt.tight_layout()
plt.savefig("espectro_PHA_asignado.png", dpi=600)
plt.show()

# === Guardar datos procesados como CSV ===
df = pd.DataFrame({'ppm': ppm_axis, 'intensidad': spectrum_smoothed})
df.to_csv("espectro_PHA_procesado.csv", index=False)
print("Espectro guardado como 'espectro_PHA_asignado.png' y 'espectro_PHA_procesado.csv'")

#########################################################################################################
# Integracion y estequiometria

# === Corrección de línea base + suavizado ===
baseline = savgol_filter(spectrum, 101, 3)
spectrum_corrected = np.clip(spectrum - baseline, 0, None)
spectrum_smoothed = savgol_filter(spectrum_corrected, 31, 2)
#print(f"la judas")
# === Eje de desplazamiento químico (ppm) ===
n_points = spectrum.shape[0]
hz_axis = np.linspace(-sw / 2, sw / 2, n_points)
ppm_axis = (hz_axis + o1) / bf1
ppm_axis = ppm_axis[::-1]
spectrum_smoothed = spectrum_smoothed[::-1]

#####8888
# Invertir para que ppm vaya de mayor a menor (como en espectros RMN)
ppm_axis = ppm_axis[::-1]
spectrum_smoothed = spectrum_smoothed[::-1]

# === RANGOS DE INTEGRACIÓN (ppm) ===
rangos = {
    "CH3_3HB": (1.15, 1.25),    # señal de CH3 (3HB)
    "CH2_3HB": (2.45, 2.60),    # señal de CH2 (3HB)
    "CH_3HB": (5.15, 5.30),     # señal de CH (3HB)
    "CH3_3HV": (0.85, 0.95),    # opcional, si hay 3HV
    "Terminal_CH3": (0.75, 0.85),  # señal terminal (opcional, para DPn)
}
#print(f"la judas")
# === Función para integrar cada región ===


# Rango valido?
#print(f"Rango ppm válido: {ppm_axis.min():.2f} a {ppm_axis.max():.2f}")

# Puntos y resolucion
print(f"Número de puntos: {len(ppm_axis)}")
print(f"Resolución aproximada (ppm): {abs(ppm_axis[1] - ppm_axis[0]):.4f}")


# Empata las longitudes


def integrar_region(ppm_axis, spectrum, rango):
    ppm_max = max(rango)
    ppm_min = min(rango)

    # Asegurar compatibilidad con eje ppm descendente
    idx = np.where((ppm_axis <= ppm_max) & (ppm_axis >= ppm_min))[0]

    if idx.size == 0:
        print(f"⚠️  Rango {rango} no contiene datos en el espectro.")
        return 0.0

    area = trapezoid(spectrum[idx], ppm_axis[idx])
    return area

# Ver las longitudes de los ejes
#print(f"Longitud ppm_axis: {len(ppm_axis)}")
#print(f"Longitud spectrum_smoothed: {len(spectrum_smoothed)}")

#Integración ===
areas = {}
for nombre, rango in rangos.items():
    areas[nombre] = integrar_region(ppm_axis, spectrum_smoothed, rango)

# === Relaciones mol/mol ===
# Normaliza con base en CH (1H) como referencia
nH = {
    "CH3_3HB": 3,
    "CH2_3HB": 2,
    "CH_3HB": 1,
    "CH3_3HV": 3,
    "Terminal_CH3": 3,
}

mol_equiv = {}
ref_area = areas["CH_3HB"] / nH["CH_3HB"]

for señal, area in areas.items():
    mol_equiv[señal] = area / nH.get(señal, 1) / ref_area

# === Estimación del grado de polimerización (DPn) ===
if "Terminal_CH3" in areas and areas["Terminal_CH3"] > 0:
    DPn = (areas["CH3_3HB"] + areas.get("CH3_3HV", 0)) / areas["Terminal_CH3"]
else:
    DPn = None

# === Resultados ===
print("\n📊 Resultados de integración (áreas aproximadas):")
for k, v in areas.items():
    print(f" - {k}: {v:.2f}")

print("\n🧪 Relaciones mol/mol normalizadas a CH (1H):")
for k, v in mol_equiv.items():
    print(f" - {k}: {v:.2f} mol")

if DPn:
    print(f"\n🔗 Estimación del grado de polimerización (DPn): {DPn:.2f}")
else:
    print("\n🔗 No se detectó señal terminal clara. No se puede estimar el DPn.")

# === Gráfica con rangos de integración ===
plt.figure(figsize=(12,6))
plt.plot(ppm_axis, spectrum_smoothed, color='black')
plt.xlim(6, 0.5)
#plt.ylim(0, max(spectrum_smoothed) * 1.2)  # Puedes ajustar el factor
plt.title("¹H-RMN de PHA: Integración por región")
plt.xlabel("ppm")
plt.ylabel("Intensidad")
colors = ["red", "blue", "green", "orange", "purple"]
for i, (label, rango) in enumerate(rangos.items()):
    plt.axvspan(rango[0], rango[1], color=colors[i%len(colors)], alpha=0.3, label=label)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("integracion_PHA.png", dpi=300)
plt.show()
