import numpy as np
import matplotlib.pyplot as plt

# Parámetros del problema
Lx, Ly = 0.1, 0.01  # Tamaño del dominio (m) (10 mm x 1 mm)
Nx, Ny = 50, 10     # Número de divisiones espaciales (aumenta Nx y Ny para mayor precisión)
dx, dy = Lx / Nx, Ly / Ny  # Tamaño de las celdas
dt = 0.001          # Paso de tiempo (s)
t_max = 10          # Tiempo máximo de simulación (s)

# Propiedades físicas
C_hielo = 2100      # Capacidad calorífica específica del hielo (J/(kg·K))
rho_hielo = 917     # Densidad del hielo (kg/m³)
k_hielo = 2.2       # Conductividad térmica del hielo (W/(m·K))
L_f = 334000        # Calor latente de fusión (J/kg)

# Inicialización de las matrices
T = np.full((Nx, Ny), -10.0)  # Temperatura inicial (°C)
phi = np.zeros((Nx, Ny))      # Fracción de fase (0: sólido, 1: líquido)

# Condiciones iniciales
T[:, -1] = np.linspace(-10, 85, Nx)  # Rampa lineal en el borde derecho

# Coeficientes de estabilidad
alpha = k_hielo / (rho_hielo * C_hielo)  # Difusividad térmica
dt_max = min(dx, dy)**2 / (4 * alpha)   # Paso de tiempo crítico (estabilidad)
if dt > dt_max:
    raise ValueError(f"dt ({dt}) supera el límite de estabilidad ({dt_max:.4e}). Reducir dt.")

# Función para calcular el siguiente paso temporal
def actualizar(T, phi):
    T_new = T.copy()
    phi_new = phi.copy()

    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            if 0 <= phi[i, j] < 1:  # Cambio de fase
                # Cambio de fracción de fase basado en la temperatura
                dT = T[i, j] - 0
                delta_phi = (dT * C_hielo) / L_f
                phi_new[i, j] = min(1, max(0, phi[i, j] + delta_phi))
                
                if phi_new[i, j] < 1:  # Mantener la temperatura en 0 °C durante el cambio de fase
                    T_new[i, j] = 0
            else:
                # Conductividad térmica (esquema explícito)
                d2T_dx2 = (T[i+1, j] - 2 * T[i, j] + T[i-1, j]) / dx**2
                d2T_dy2 = (T[i, j+1] - 2 * T[i, j] + T[i, j-1]) / dy**2
                T_new[i, j] += alpha * dt * (d2T_dx2 + d2T_dy2)
    
    return T_new, phi_new

# Simulación
tiempos = np.arange(0, t_max, dt)
resultados_T = []
resultados_phi = []

for t in tiempos:
    T, phi = actualizar(T, phi)

    # Guardar resultados cada segundo
    if int(t * 1000) % 1000 == 0:
        resultados_T.append(T.copy())
        resultados_phi.append(phi.copy())

# Visualización
fig, axs = plt.subplots(2, len(resultados_T), figsize=(20, 8))
for idx, (T_plot, phi_plot) in enumerate(zip(resultados_T, resultados_phi)):
    # Gráfica de temperatura
    im1 = axs[0, idx].imshow(T_plot, cmap="coolwarm", origin="lower", extent=[0, Lx * 100, 0, Ly * 100])
    axs[0, idx].set_title(f"Temperatura t={idx}s")
    axs[0, idx].set_xlabel("Posición (cm)")
    axs[0, idx].set_ylabel("Altura (cm)")
    fig.colorbar(im1, ax=axs[0, idx])

    # Gráfica de fracción de fase
    im2 = axs[1, idx].imshow(phi_plot, cmap="viridis", origin="lower", extent=[0, Lx * 100, 0, Ly * 100])
    axs[1, idx].set_title(f"Fracción de fase t={idx}s")
    axs[1, idx].set_xlabel("Posición (cm)")
    axs[1, idx].set_ylabel("Altura (cm)")
    fig.colorbar(im2, ax=axs[1, idx])

plt.tight_layout()
plt.show()