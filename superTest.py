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

# Definir las regiones de agua (horizontal) y hielo (vertical)
""" T[:, :Ny//4] = T_agua  # Agua en la parte izquierda
T[:, 3*Ny//4:] = T_agua  # Agua en la parte derecha """

# Establecer los 10 nodos centrales a -20
centro_inicio = Nx//2 - 5
centro_fin = Nx//2 + 5

#Construccion de Matriz de temperaturas
T = np.ones((Nx, Nx)) * T_inicial  # Todo el dominio como hielo inicialmente
T[:, 0:centro_inicio] = 20  
T[:, centro_fin:] = 20
T[Ny:,:] = np.nan
T[Ny:, centro_inicio:centro_fin] = 20

#Construccion de Matriz de fraccion de fase
phi = np.zeros((Nx, Nx))  # Fracción de fase inicial
phi[:, 0:centro_inicio] = 1
phi[:, centro_fin:] = 1
phi[Ny:, :] = np.nan
phi[Ny:, centro_inicio:centro_fin] = 1

# Precomputar constantes
alpha = k / (rho * c)
dt = 0.25 * min(dx**2, dy**2) / alpha  # Ajuste según estabilidad explícita

factor_x = alpha * dt / dx**2
factor_y = alpha * dt / dy**2 

# Verificar estabilidad numérica
Fo_x = factor_x * dx**2 / alpha
Fo_y = factor_y * dy**2 / alpha
if Fo_x + Fo_y > 0.5:
    raise ValueError("El paso de tiempo dt no cumple el criterio de estabilidad.")

# Función para actualizar fracción de fase
def actualizar_phi(T, phi, L, c):
    for j in range(phi.shape[0]):
        for i in range(phi.shape[1]):
            if 0 <= phi[j, i] < 1 and phi[j, i] != None:  # 0 todo hielo, 1 todo agua
                dT = T[j, i] - T_m
                d_phi = (dT * c) / L
                phi[j, i] = min(1, max(0, phi[j, i] + d_phi))  # Mantiene phi en rango
                if phi[j, i] < 1 and T[j, i] >= T_m:
                    T[j, i] = T_m
    return phi

# Simulación
tiempo_total = 0
tiempos = []
phi_promedios = []
T_sol = [T]

t = 0.0
t_final = 1

while t <= t_final:
    T_dt = np.copy(T)
    
    # Aplicar diferencias finitas explícitas solo en la parte superior de la "T"
    d2T_dx2 = factor_x * (np.roll(T[:Ny, centro_inicio:centro_fin], -1, axis=1) - 2 * T[:Ny, centro_inicio:centro_fin] + np.roll(T[:Ny, centro_inicio:centro_fin], 1, axis=1)) / dx**2
    d2T_dy2 = factor_y * (np.roll(T[:Ny, centro_inicio:centro_fin], -1, axis=0) - 2 * T[:Ny, centro_inicio:centro_fin] + np.roll(T[:Ny, centro_inicio:centro_fin], 1, axis=0)) / dy**2 

    T_dt[:Ny, centro_inicio:centro_fin] += dt * (d2T_dx2 + d2T_dy2)

    #print(d2T_dx2)

    # Condiciones de borde (aislamiento térmico → derivada nula)
    T_dt[:Ny, centro_inicio-1] = 20
    T_dt[:Ny, centro_fin+1] = 20
    
    """ T_dt[0, :] = T_dt[1, :]     # Borde superior
    T_dt[Ny-1, :] = T_dt[Ny-2, :]  # Borde inferior de la parte superior
    T_dt[:, 0] = T_dt[:, 1]     # Borde izquierdo
    T_dt[:, -1] = T_dt[:, -2]   # Borde derecho """
    
    T = np.copy(T_dt)
    
    """plt.imshow(T, cmap='coolwarm', origin='lower', extent=[0, Lx*1e3, 0, Lx*1e3])
    plt.colorbar(label="Temperatura (°C)")
    plt.xlabel("Longitud (mm)")
    plt.ylabel("Altura (mm)")
    plt.title("Distribución de temperatura en el dominio")
    plt.show()"""
    
    t += dt

    

#print(T[:Ny, :])
# Graficar la temperatura final
plt.imshow(phi, cmap='coolwarm', origin='lower', extent=[0, Lx*1e3, 0, Lx*1e3])
plt.colorbar(label="Temperatura (°C)")
plt.xlabel("Longitud (mm)")
plt.ylabel("Altura (mm)")
plt.title("Distribución de temperatura en el dominio")
plt.show()




""" while True:
    # Actualizar fracción de fase
    phi = actualizar_phi(T)

    # Actualizar temperatura (esquema explícito vectorizado)
    T[1:-1, 1:-1] += factor_x * (T[2:, 1:-1] - 2 * T[1:-1, 1:-1] + T[:-2, 1:-1]) + \
                     factor_y * (T[1:-1, 2:] - 2 * T[1:-1, 1:-1] + T[1:-1, :-2])

    # Condiciones de borde (aislamiento térmico)
    T[0, :] = T[1, :]       # Borde superior
    T[-1, :] = T[-2, :]     # Borde inferior
    T[:, 0] = T[:, 1]       # Borde izquierdo
    T[:, -1] = T[:, -2]     # Borde derecho

    tiempo_total += dt
    T_sol.append(T.copy())
    
    # Verificar criterio de parada en la sección vertical
    phi_promedio = np.mean(phi[Nx//4:3*Nx//4, Ny//2])  # Promedio en la sección vertical
    tiempos.append(tiempo_total)
    phi_promedios.append(phi_promedio)

    if phi_promedio >= 0.5:
        break

    # Visualización periódica (opcional)
    if int(tiempo_total) % 1 == 0:  # Graficar cada segundo
        plt.imshow(T, cmap='coolwarm', origin='lower')
        plt.colorbar(label='Temperatura (°C)')
        plt.title(f'Tiempo: {tiempo_total:.2f} s')
        plt.pause(0.1)

T_array = np.array(T_sol)
# Mostrar el tiempo final
print(f"Tiempo para recuperar 50% de sección útil: {tiempo_total:.2f} s") """
# Graficar Temperaturas
""" im = plt.imshow(T_array.T, aspect='auto', cmap='hot', origin='lower',
                 extent=[0, tiempo_total, 0, Ly])

# Graficar la evolución de la fracción de fase promedio
plt.figure()
plt.plot(tiempos, phi_promedios, label="Fracción de fase promedio")
plt.axhline(0.5, color='r', linestyle='--', label="50% de fracción de fase")
plt.xlabel("Tiempo (s)")
plt.ylabel("Fracción de fase promedio")
plt.title("Evolución de la fracción de fase en la sección vertical")
plt.legend()
plt.grid()
plt.show() """
