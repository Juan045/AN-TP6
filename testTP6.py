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
T = np.ones(cantidad_divisiones) #Se utiliza como vector de temperaturas de la barra

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
t_sol = [t]

#METODO IMPLICITO
while t < t_final:
    A = np.zeros((cantidad_divisiones, cantidad_divisiones))
    b = np.zeros(cantidad_divisiones)
    for i in range(cantidad_divisiones):
        if i == 0:
            A[i][i] = 1
            b[i] = T0
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
    T_dt = np.dot(np.linalg.inv(A), b)
    
    T = T_dt.copy()
    t += dt
    T_sol.append(T)
    t_sol.append(t)
    
    if int(t * 1000) % int(10 * 1000) == 0:
        plt.plot(x, T, label=f"t = {int(t)} s")
        
plt.xlabel("Posición (mm)")
plt.ylabel("Temperatura (°C)")
plt.title("Evolución de la temperatura en el dominio")
plt.legend()
plt.grid(True)
plt.show()  
print(t)

    
# GRAFICAS
""" import matplotlib.animation as animation
fig = plt.figure()
ax = plt.gca()

def actualizar(i):
    ax.clear()
    plt.plot(x, T_sol[i], "ro")
    plt.title(str(round(t_sol[i], 5)))
    plt.xlim(0, L)
    plt.ylim(-20, 100)
    
animacion = animation.FuncAnimation(fig, actualizar, range(len(t_sol)))
plt.show() """

