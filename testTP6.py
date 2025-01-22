import numpy as np
import matplotlib.pyplot as plt

def rampa_lineal(t):
    if (t <= 10):
        return -10 + 9.5*t
    else:
        return 85.

#DATOS
#Se desprecia el ancho de 1 mm del dominio
L = 1.0 #mm largo del dominio
C_agua = 1 #cal/(g*°C)
ro_agua = 0.998 #g/(cm^3) -> es 1 a 4°C
C_hielo = 0.5 #cal/(g*°C)
ro_hielo = 0.917 #g/(cm^3)
k_prima_agua = 0.0014 #cal/(s*cm*°C)
k_prima_hielo = 0.005 #cal/(s*cm*°C)
L_f = 80  # cal/g (calor latente de fusión)

k_agua = k_prima_agua/(ro_agua*C_agua)
k_hielo = k_prima_hielo/(ro_hielo*C_hielo)

print(k_agua, k_hielo)

dx = 0.1 #Test con 1mm
dt = 0.5 * (dx ** 2) #Establecer un criterio para justificarlo (cambiar si es necesario)
dt = 0.001

landa_agua = k_agua * (dt/dx**2)
landa_hielo = k_hielo * (dt/dx**2)

#DISCRETIZACION DEL DOMINIO
cantidad_divisiones = int(L/dx)
x = np.linspace(0, L, cantidad_divisiones)

#CONDICIONES INICIALES
t = 0
t_inicial = t
T = np.ones(cantidad_divisiones) * -10 #Se utiliza como vector de temperaturas de la barra
phi = np.zeros(cantidad_divisiones)  #Fracción de fase inicial (todo es hielo)

#CONDICIONES DE BORDE
T0 = -10 #como esta aislado considero 0, temperatura el inicio del dominio
T_F = rampa_lineal(0) #Temperatura al final del dominio
T[0] = T0
T[-1] = T_F

print(f"Delta estable: {(dx/10)**2/(2*k_hielo)} -> dt = {dt}")
print(f"Delta estable: {(dx/10)**2/(2*k_agua)}")

#SOLUCIONES DEL SISTEMA
t_final = 100 #seg
T_dt = T.copy()
T_sol = [T]
phi_sol = [phi]
t_sol = [t]

#METODO IMPLICITO
while t < 100:
    A = np.zeros((cantidad_divisiones, cantidad_divisiones))
    b = np.zeros(cantidad_divisiones)
    for i in range(cantidad_divisiones):
        if i == 0:
            A[i][i] = 1
            b[i] = T[i + 1]
        elif i == cantidad_divisiones - 1:
            A[i][i] = 1
            b[i] = rampa_lineal(t)
        else:
            if T[i] < 0:
                lambda_ = landa_hielo
            else:
                lambda_ = landa_agua
            A[i][i] = 1 + 2*lambda_
            A[i][i-1]=-lambda_
            A[i][i+1]=-lambda_
            b[i] = T[i]
            
    T_dt = np.dot(np.linalg.inv(A), b) #Resolver nuevas temperaturas
    
    #Corregir temperatura
    for i in range(cantidad_divisiones):
        if 0 <= phi[i] < 1:
            dT = T_dt[i] - T[i]
            d_phi = (dT * C_hielo) / L_f
            phi[i] = min(1, max(0, phi[i] + d_phi)) #Mantiene el valor de phi dentro del rango de 0 a 1
            if phi[i] < 1:
                T_dt[i] = 0
    
    t += dt
    T = T_dt.copy() #Nuevas temperaturas
    
    if int(t * 1000) % int(10 * 1000) == 0:
        T_sol.append(T)
        phi_sol.append(phi.copy())
        t_sol.append(t)

#MOSTRAR RESULTADOS FINALES
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

#Gráfico de temperatura
for i, t in enumerate(t_sol):
    ax1.plot(x, T_sol[i], label=f"T = {int(t)} s")
ax1.set_title("Evolución de la temperatura")
ax1.set_xlabel("Posición (cm)")
ax1.set_ylabel("Temperatura (°C)")
ax1.legend()
ax1.grid(True)

#Gráfico de fracción de fase
for i, t in enumerate(t_sol):
    ax2.plot(x, phi_sol[i], label=f"Fase t = {int(t)} s")
ax2.set_title("Evolución de la fracción de fase")
ax2.set_xlabel("Posición (cm)")
ax2.set_ylabel("Fracción de fase")
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()
    
    
""" if int(t * 1000) % int(10 * 1000) == 0:
    plt.plot(x, T, label=f"t = {int(t)} s")
    plt.plot(x, phi, label=f"Fase t = {int(t)} s") """
        
""" plt.xlabel("Posición (mm)")
plt.ylabel("Temperatura (°C)")
plt.title("Evolución de la temperatura en el dominio")
plt.legend()
plt.grid(True)
plt.show()  
print(t) """

    
""" # Animación de los gráficos
from matplotlib.animation import FuncAnimation
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

def actualizar(frame):
    ax1.clear()
    ax2.clear()

    # Gráfico de temperatura
    ax1.plot(x, T_sol[frame], color="red")
    ax1.set_title("Evolución de la temperatura")
    ax1.set_xlabel("Posición (cm)")
    ax1.set_ylabel("Temperatura (°C)")
    ax1.grid(True)

    # Gráfico de fracción de fase
    ax2.plot(x, phi_sol[frame], color="blue")
    ax2.set_title("Evolución de la fracción de fase")
    ax2.set_xlabel("Posición (cm)")
    ax2.set_ylabel("Fracción de fase")
    ax2.grid(True)

    fig.tight_layout()

# Crear la animación
anim = FuncAnimation(fig, actualizar, frames=len(t_sol), interval=50)
plt.show() """

