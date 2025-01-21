import numpy as np
import explicito
import matplotlib.pyplot as plt

# DATOS
L=10 #Longitud de la barra (cm)
Kprima=0.49
C=0.2174
ro=2.7
k=Kprima/(ro*C)
deltax=2 #(cm)
dt=0.1 #(seg)
landa=k*(dt/deltax**2)

# DISCRETIZACION DE LA BARRA
N=int(L/deltax) #Cuantas divisiones se aplicaran a la barra
x=np.linspace(0,L,N) #Las divisiones de la barra

# CONDICION INICIAL
t=0
Tinicial=0
T=np.ones(N) #Vector de temperaturas de la barra

# CONDICIONES DE BORDE
T0=100 #Temperatura al inicio de la barra
TL=50 #Temperatura al final de la barra
T[0]=T0
T[-1]=TL

# CONVERGENCIA Y ESTABILIDAD
deltaTestable=deltax**2/2*k # Para que el sistema sea estable y converga debe cumplirse esta condicion
# deltaT entonces => 1.66995 aprox.

# SOLUCIONES DEL SISTEMA
tfinal=12 # (seg)
Tdt=T.copy() #Temperatura del siguiente instante de tiempo
Tsol=[T] #Lista donde cada elemento es un paso de tiempo que guarda los vectores de temperatura
tsol=[t]

# CALCULO DE LA TEMPERATURA
# Metodo explicito
""" while t<tfinal:
    for i in range(N):
        if i==0:
            Tdt[i]=T0 #Condicion de borde
        elif i==N-1:
            Tdt[i]=TL #Condicion de borde
        else:
            Tdt[i]=T[i]+landa*(T[i+1]-(2*T[i])+(T[i-1])) #Formula de ec. diferencial finita
    T=Tdt.copy() #Actualiza el vector de temperaturas para la proxima iteracion
    t=t+dt
    Tsol.append(T)
    tsol.append(t) """
    
#Metodo Implicito
while t<tfinal:
    A=np.zeros((N,N))
    b=np.zeros(N)
    for i in range(N):
        if i==0:
            A[i][i]=1 #Primer elemento de la matriz
            b[i]=T0
        elif i==N-1:
            A[i][i]=1 #Ultimo elemento de la matriz
            b[i]=TL
        else: #ec. diferencial finita
            A[i][i]=1+2*landa
            A[i][i-1]=-landa
            A[i][i+1]=-landa
            b[i]=T[i]
    #Resolver el sistema de ecuaciones
    Ainv=np.linalg.inv(A)
    Tdt=np.dot(Ainv,b)
    
    T=Tdt.copy()
    t=t+dt
    Tsol.append(T)
    tsol.append(t)


# GRAFICAS
i=10
# ANIMACION       
import matplotlib.animation as animation
fig=plt.figure()
ax=plt.gca()

def actualizar(i):
    ax.clear()
    plt.plot(x,Tsol[i],"ro")
    plt.title(str(round(tsol[i],5)))
    plt.xlim(0,L)
    plt.ylim(0,100)
ani=animation.FuncAnimation(fig,actualizar,range(len(tsol)))
plt.show()
