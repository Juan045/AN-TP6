import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos
dx = 0.1  # Tamaño de celda [mm]
dt = 0.01  # Paso de tiempo [s]
L = 10  # Longitud del hielo [mm]
W = 1  # Ancho del agua [mm]

k_hielo = 2.2  # Conductividad térmica del hielo [W/(m·K)]
k_agua = 0.6  # Conductividad térmica del agua [W/(m·K)]
C_hielo = 2000  # Capacidad calorífica del hielo [J/(kg·K)]
C_agua = 4186  # Capacidad calorífica del agua [J/(kg·K)]
rho_hielo = 917  # Densidad del hielo [kg/m³]
rho_agua = 1000  # Densidad del agua [kg/m³]
L_f = 334000  # Calor latente de fusión [J/kg]

# Dimensiones del dominio
nx_hielo = int(L / dx)
ny_hielo = int(W / dx)
nx_total = nx_hielo
ny_total = 2 * ny_hielo + nx_hielo

dominio = np.zeros((nx_total, ny_total))  # Dominio 2D
T = np.full((nx_total, ny_total), -10.0)  # Temperatura inicial [°C]
phi = np.zeros((nx_total, ny_total))  # Fracción de fase inicial

# Inicialización del dominio
# Región de hielo
for i in range(nx_hielo):
    for j in range(ny_hielo, ny_hielo + nx_hielo):
        phi[i, j] = 0  # Hielo

# Región de agua
for i in range(nx_hielo):
    for j in range(ny_total):
        if j < ny_hielo or j >= ny_hielo + nx_hielo:
            phi[i, j] = 1  # Agua
            T[i, j] = 0

# Función para calcular el coeficiente de difusión

def alpha(phi, k_hielo, k_agua, rho_hielo, rho_agua, C_hielo, C_agua):
    k = phi * k_agua + (1 - phi) * k_hielo
    rho = phi * rho_agua + (1 - phi) * rho_hielo
    C = phi * C_agua + (1 - phi) * C_hielo
    return k / (rho * C)

# Bucle de simulación
num_steps = 1000  # Número de pasos de tiempo
for step in range(num_steps):
    T_new = T.copy()
    for i in range(1, nx_total - 1):
        for j in range(1, ny_total - 1):
            # Coeficiente de difusión
            a = alpha(phi[i, j], k_hielo, k_agua, rho_hielo, rho_agua, C_hielo, C_agua)

            # Ecuación de calor (diferencias finitas)
            d2T_dx2 = (T[i + 1, j] - 2 * T[i, j] + T[i - 1, j]) / dx**2
            d2T_dy2 = (T[i, j + 1] - 2 * T[i, j] + T[i, j - 1]) / dx**2
            T_new[i, j] = T[i, j] + dt * a * (d2T_dx2 + d2T_dy2)

            # Cambio de fase
            if 0 <= phi[i, j] < 1:
                dT = T_new[i, j] - T[i, j]
                d_phi = (dT * C_hielo * rho_hielo) / L_f
                phi[i, j] = min(1, max(0, phi[i, j] + d_phi))
                if phi[i, j] < 1 and T_new[i, j] >= 0:
                    T_new[i, j] = 0

    T = T_new.copy()

    # Visualización en cada 50 pasos
    if step % 50 == 0:
        plt.figure(figsize=(8, 6))
        plt.imshow(phi, cmap="coolwarm", origin="lower", extent=[0, ny_total * dx, 0, nx_total * dx])
        plt.colorbar(label="Fracción de agua (phi)")
        plt.title(f"Paso de tiempo {step}")
        plt.xlabel("mm (ancho)")
        plt.ylabel("mm (largo)")
        plt.show()
