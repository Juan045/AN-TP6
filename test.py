import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# Parámetros
Lx, Ly = 1.0, 1.0  # Dimensiones del dominio (cm)
dx, dy = 0.05, 0.05  # Tamaño de celdas (cm)
nx, ny = int(Lx / dx), int(Ly / dy)  # Número de divisiones
dt = 0.001  # Paso de tiempo (s)
t_final = 100  # Tiempo de simulación (s)

C_hielo, C_agua = 0.5, 1.0  # Capacidad calorífica específica (cal/g°C)
rho_hielo, rho_agua = 0.917, 0.998  # Densidad (g/cm³)
k_hielo, k_agua = 0.005 / (rho_hielo * C_hielo), 0.0014 / (rho_agua * C_agua)
L_f = 79.8279  # Calor latente (cal/g)

# Difusividad térmica
alpha_hielo = k_hielo
alpha_agua = k_agua

# Inicialización del dominio
T = np.ones((nx, ny)) * -10  # Temperatura inicial (°C)
phi = np.zeros((nx, ny))  # Fracción de fase inicial (hielo)

# Bordes
T[:, -1] = np.linspace(-10, 85, nx)  # Rampa de temperatura en el borde derecho

# Matriz de coeficientes para el método implícito
def construir_matriz(nx, ny, alpha, dx, dy, dt):
    N = nx * ny
    diagonales = []
    posiciones = []

    # Diagonal principal
    diagonales.append((1 + 2 * alpha * dt / dx**2 + 2 * alpha * dt / dy**2) * np.ones(N))
    posiciones.append(0)

    # Diagonales para las interacciones en x
    diagonales.append((-alpha * dt / dx**2) * np.ones(N - 1))
    diagonales.append((-alpha * dt / dx**2) * np.ones(N - 1))
    posiciones.extend([-1, 1])

    # Diagonales para las interacciones en y
    diagonales.append((-alpha * dt / dy**2) * np.ones(N - nx))
    diagonales.append((-alpha * dt / dy**2) * np.ones(N - nx))
    posiciones.extend([-nx, nx])

    return diags(diagonales, posiciones, shape=(N, N))

# Aplanar la matriz T para trabajar con el método implícito
T_flat = T.flatten()

# Construir la matriz A
alpha = alpha_hielo  # Usamos hielo por defecto
A = construir_matriz(nx, ny, alpha, dx, dy, dt).tocsc()

# Simulación
t = 0
while t < t_final:
    if int(t * 1000) % int(10 * 1000) == 0:
        print(t)
    # Resolver el sistema lineal
    T_new_flat = spsolve(A, T_flat)
    T_new = T_new_flat.reshape((nx, ny))

    # Actualizar fracción de fase
    for i in range(nx):
        for j in range(ny):
            if 0 <= phi[i, j] < 1:
                dT = T_new[i, j] - T[i, j]
                d_phi = (dT * C_hielo) / L_f
                phi[i, j] = min(1, max(0, phi[i, j] + d_phi))
                if phi[i, j] < 1 and T_new[i, j] >= 0:
                    T_new[i, j] = 0

    # Actualizar temperatura y tiempo
    T = T_new.copy()
    T_flat = T.flatten()
    t += dt

# Visualización final
plt.imshow(T, cmap='hot', origin='lower', extent=[0, Lx, 0, Ly])
plt.colorbar(label="Temperatura (°C)")
plt.title("Distribución de Temperatura al Final")
plt.xlabel("x (cm)")
plt.ylabel("y (cm)")
plt.show()