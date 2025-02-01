import numpy as np
import matplotlib.pyplot as plt

import sys

#DATOS 
C_agua = 1 #cal/(g*°C)
ro_agua = 0.998 #g/(cm^3) -> es 1 a 4°C
C_hielo = 0.5 #cal/(g*°C)
ro_hielo = 0.917 #g/(cm^3)
k_prima_agua = 0.0014 #cal/(s*cm*°C)
k_prima_hielo = 0.005 #cal/(s*cm*°C)
L_f = 79.8279  # cal/g (calor latente de fusión)

k_agua = k_prima_agua/(ro_agua*C_agua)
k_hielo = k_prima_hielo/(ro_hielo*C_hielo)

#DOMINIO
Lx = 10e-3
Ly = 1e-3 
Nx = 30   
Ny = 10
dx = Lx / Nx
dy = Ly / Ny

#CONDICIONES INICIALES
T_hielo = -20
T_agua = 20

# Establecer los 10 nodos centrales a -20
centro_inicio = Nx//2 - 5
centro_fin = Nx//2 + 5

#Construccion de Matriz de temperaturas
T = np.ones((Nx, Nx)) * T_hielo  # Todo el dominio como hielo inicialmente
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

print(f"El paso de tiempo recomendado es: {0.25 * min(dx** 2, dy** 2) / k_hielo}")
dt = 0.5 * min(dx**2, dy**2) / min(k_hielo, k_agua)
landa_agua = k_agua * (dt/dx**2)
landa_hielo = k_hielo * (dt/dx**2)

t = 0.0
t_final = 0.001


print(f"dt = {dt}")


while t < t_final:
    T_dt = np.copy(T)
    
    T_dx2 = (np.roll(T_dt, -1, axis=1) - 2 * T_dt + np.roll(T_dt, 1, axis=1)) / dx**2
    #T_dy2 = (np.roll(T_dt, -1, axis=0) - 2 * T_dt + np.roll(T_dt, 1, axis=0)) / dy**2
                 
    #T_dt[:Ny, centro_inicio+1:centro_fin] += landa_hielo * (T_dt[:Ny, centro_inicio:centro_fin-2] - 2 * T_dt[:Ny, centro_inicio+1:centro_fin] + T_dt[:Ny, centro_inicio+2:centro_fin])
    #T_dt[1:Ny, centro_inicio:centro_fin] += landa_hielo * (T_dt[:Ny-2, centro_inicio:centro_fin] - 2 * T_dt[1:Ny, centro_inicio:centro_fin] + T_dt[2:Ny, centro_inicio:centro_fin])             
                    
    T_dt += landa_hielo * dt*0.001 * (T_dx2)
    
    #CONDICIONES DE BORDE
    #Rama superior
    T_dt[0,:] = T_dt[1, :]
    T_dt[:,0] = T_dt[:,1]
    T_dt[:,-1] = T_dt[:,-2]
    #Rama superior izq
    T_dt[Ny, :centro_inicio] = T_dt[Ny - 1, :centro_inicio] 
    #Rama superior der
    T_dt[Ny, centro_fin:] = T_dt[Ny - 1, centro_fin:] 
    #Rama inferior, lado izquierdo
    T_dt[Ny:, centro_inicio] = T_dt[Ny:, centro_inicio+1]
    #Rama inferior, lado derecho
    T_dt[Ny:, centro_fin-1] = T_dt[Ny:, centro_fin-2] 
    
    T = np.copy(T_dt) 
    t += dt



plt.imshow(T, cmap='coolwarm', origin='lower', extent=[0, Lx*1e3, 0, Lx*1e3])
plt.colorbar(label="Temperatura (°C)")
plt.xlabel("Longitud (mm)")
plt.ylabel("Altura (mm)")
plt.title("Distribución de temperatura en el dominio")
plt.show()


plt.imshow(phi, cmap='coolwarm', origin='lower', extent=[0, Lx*1e3, 0, Lx*1e3])
plt.colorbar(label="Fraccion")
plt.xlabel("Longitud (mm)")
plt.ylabel("Altura (mm)")
plt.title("Distribución de fracción de fase en el dominio")
plt.show()