import numpy as np
import matplotlib.pyplot as plt

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

#CANALES DE AGUA
L_a = 1.0 #(cm)
##Condiciones de borde
T_0a = 20 #[°C]
T_Fa = 20 #[°C]

#TROZO DE HIELO
L_h = 0.1 #[cm]
W_h = 0.1 #[cm]
##Condiciones de borde
T_0h = -10 #[°C]
T_Fh = -10 #[°C]

#DISCRETIZACIONES
dx = 0.05
dy = 0.05
dt = 0.001

##VER ESTO
landa_agua = k_agua * (dt/dx**2)
landa_hielo = k_hielo * (dt/dx**2) 