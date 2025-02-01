import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos
rho = 917  # kg/m³ (densidad del hielo)
c = 2100   # J/(kg·K) (calor específico)
k = 2.2    # W/(m·K) (conductividad térmica)
L = 334000 # J/kg (calor latente de fusión)
T_m = 0    # °C (temperatura de fusión)
delta_T = 0.1  # Rango de transición de fase

# Parámetros del dominio
Lx, Ly = 11e-3, 11e-3  # Dimensiones en metros (11x11 mm)
Nx, Ny = 110, 110      # Número de nodos (ajustar según resolución deseada)
dx, dy = Lx / Nx, Ly / Ny

# Condiciones iniciales
T_inicial = -20  # °C (temperatura inicial del hielo)
T_agua = 20      # °C (temperatura del agua)

# Inicialización del campo de temperatura
T = np.ones((Ny, Nx)) * T_inicial  # Todo el dominio como hielo inicialmente
phi = np.zeros((Ny, Nx))  # Fracción de fase inicial

# Definir regiones de agua (T-shape)
# Agua en el fondo (10 mm)
T[-int(10e-3 / dy):, :] = T_agua
phi[-int(10e-3 / dy):, :] = 1

# Agua en los lados (5 mm a cada lado)
T[:, :int(5e-3 / dx)] = T_agua
T[:, -int(5e-3 / dx):] = T_agua
phi[:, :int(5e-3 / dx)] = 1
phi[:, -int(5e-3 / dx):] = 1

# Definir el bloque de hielo (1x1 mm) en la parte superior central
ice_start_x = int((Nx - int(1e-3 / dx)) / 2)
ice_end_x = ice_start_x + int(1e-3 / dx)
ice_start_y = int((Ny - int(1e-3 / dy)) / 2) - int(5e-3 / dy)  # Centrado en la parte superior
ice_end_y = ice_start_y + int(1e-3 / dy)

T[ice_start_y:ice_end_y, ice_start_x:ice_end_x] = T_inicial
phi[ice_start_y:ice_end_y, ice_start_x:ice_end_x] = 0

# Precomputar constantes
alpha = k / (rho * c)
dt = 0.25 * min(dx**2, dy**2) / alpha  # Ajuste según estabilidad explícita
factor_x = alpha * dt / dx**2
factor_y = alpha * dt / dy**2

# Función para actualizar la fracción de fase
def actualizar_phi(T, phi, L, c):
    for j in range(phi.shape[0]):
        for i in range(phi.shape[1]):
            if 0 <= phi[j, i] < 1:  # 0 todo hielo, 1 todo agua
                dT = T[j, i] - T_m
                d_phi = (dT * c) / L
                phi[j, i] = min(1, max(0, phi[j, i] + d_phi))  # Mantiene phi en rango
                if phi[j, i] < 1 and T[j, i] >= T_m:
                    T[j, i] = T_m
    return phi

# Simulación
t = 0
tiempo_max = 100  # segundos

tiempos = []
phi_promedios = []

while t < tiempo_max:
    T_new = np.copy(T)

    # Aplicar ecuación de calor explícita (Diferencias Finitas)
    T_new[:, 1:-1] = T[:, 1:-1] + factor_x * (T[:, :-2] - 2*T[:, 1:-1] + T[:, 2:])
    T_new[1:-1, :] += factor_y * (T[:-2, :] - 2*T[1:-1, :] + T[2:, :])

    # Actualizar fracción de fase
    phi = actualizar_phi(T_new, phi, L, c)
    phi_promedio = np.mean(phi[ice_start_y:ice_end_y, ice_start_x:ice_end_x])  # Promedio de la fracción de fase en el hielo

    # Actualizar temperatura
    T = np.copy(T_new)
    t += dt
    tiempos.append(t)
    phi_promedios.append(phi_promedio)

    if phi_promedio >= 0.5:
        break

# Graficar resultados
print(f"Tiempo para que el 50% del hielo se derrita: {t} segundos")
plt.figure(figsize=(8, 8))
plt.imshow(T, cmap='coolwarm', aspect='auto', origin='lower', extent=[0, Lx*1e3, 0, Ly*1e3])
plt.colorbar(label="Temperatura (°C)")
plt.xlabel("Longitud (mm)")
plt.ylabel("Altura (mm)")
plt.title("Distribución de Temperatura en el Dominio")
plt.show()
