#PROCESAMIENTO DE DATOS CRUDOS RMN - IDENTIFFICACION Y CUNTIFICACIO DE SEÑLES DE PHA
#Script creado por ANDREA TREJO y LUIS MARIO HERNANDEZ-SOTO 05/2025

#CARGAR LIBRERIAS
import nmrglue as ng
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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

