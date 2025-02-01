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
Lx, Ly = 10e-3, 1e-3  # Dimensiones en metros
Nx, Ny = 110, 10      # Número de nodos
dx, dy = Lx / Nx, Ly / Ny

# Condiciones iniciales
T_inicial = -20  # °C
T_agua = 20      # °C

# Inicialización del campo de temperatura
T = np.ones((Ny, Nx)) * T_inicial  # Todo el dominio como hielo inicialmente
phi = np.zeros((Ny, Nx))  # Fracción de fase inicial

# Definir regiones de agua e hielo
centro_inicio = Nx//2 - 5
centro_fin = Nx//2 + 5
T[:, 0:centro_inicio] = T_agua
T[:, centro_fin:] = T_agua
T[-1, :] = T_agua
phi[:, 0:centro_inicio] = 1
phi[:, centro_fin:] = 1

# Precomputar constantes
alpha = k / (rho * c)
dt = 0.25 * min(dx**2, dy**2) / alpha  # Ajuste según estabilidad explícita
factor_x = alpha * dt / dx**2
factor_y = alpha * dt / dy**2

# Función para actualizar la fracción de fase
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
print(dt)

import sys
while t < tiempo_max:
    T_new = np.copy(T)
    
    print(T[:, :-2].shape)
    print(T[:, 1:-1].shape)
    sys.exit()
    # Aplicar ecuación de calor explícita (Diferencias Finitas)
    T_new[:, 1:-1] = T[:, 1:-1] + factor_x * (T[:, :-2] - 2*T[:, 1:-1] + T[:, 2:])
    T_new[1:-1, :] += factor_y * (T[:-2, :] - 2*T[1:-1, :] + T[2:, :])
    
    # Actualizar fracción de fase
    phi = actualizar_phi(T_new, phi, L, c)
    phi_promedio = np.mean(phi[:, centro_inicio:centro_fin])
    
    # Actualizar temperatura
    T = np.copy(T_new)
    t += dt
    tiempos.append(t)
    #phi_promedios.append(np.mean(phi))
    if phi_promedio >= 0.5:
        break

# Graficar resultados
print(phi[:, centro_inicio:centro_fin])
print(T[:, centro_inicio:centro_fin])
print(t)
plt.figure(figsize=(8, 4))
plt.imshow(phi[:, centro_inicio:centro_fin], cmap='coolwarm', aspect='auto', origin='lower', extent=[0, Ly*1e3, 0, Ly*1e3])
plt.colorbar(label="Temperatura (°C)")
plt.xlabel("Longitud (mm)")
plt.ylabel("Altura (mm)")
plt.title("Distribución de Temperatura en el Dominio")
plt.show()