#Aqui ira el programa que lea el fichero de estaciones

import matplotlib.pyplot as plt
import numpy as np
import utils_tempos as ut

#Extraemos de los ficheros las posiciones de las islas y de las esatciones sismicas
islas_long,islas_lat = np.loadtxt("contorno_islas.txt", unpack=True)
est_lat,est_long,est_alturas=np.loadtxt("fichero_estaciones.txt",comments='---',usecols=(4,5,6),unpack=True)
vel,capas=np.loadtxt("canary.txt",comments='!',usecols=(3,4),unpack=True)

rango_longitudes=np.arange(-17,-16,0.05)
rango_latitudes=np.arange(28.6,27.9,-0.05)
rango_profundidades=np.arange(0,30,5)

puntos_mapa=np.array([rango_longitudes,rango_latitudes,rango_profundidades])
estaciones=np.array([est_long,est_lat,est_alturas])

"""

Primero vamos a calcular las distancias de cada punto en el rango elegido a las estaciones.
Para ordenar los datos los colocaremos en una matriz 2D en la que cada columna representara
una estacion y cada fila la distancia a dicha estacion

"""

def distancia(estaciones,puntos_mapa):
    d=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1]),len(estaciones[0])))
    for k in np.arange(len(estaciones[0])):#Este bucle recorre cada estacion
        c=0
        for i in np.arange(len(puntos_mapa[1])): #Este buble recorre cada longitud (cada fila)
            for j in np.arange(len(puntos_mapa[0])): #Este bucle recorre cada latitud (cada columna dentro de cada fila)
                #Calculamos la distancia entre la estación k-esima y el punto del mapa [i][j] 
                d[c][k]=ut.distance(estaciones[0][k],estaciones[1][k],puntos_mapa[0][j],puntos_mapa[1][i])
                c=c+1
    return d
    
"""

A continuación calcularemos los tiempos de llegada de las ondas desde cada punto del rango elegido a cada estacion.
Como en el caso de las distancias, organizamos los datos en una matriz 2D, en la que cada columna
se corresponde con una estacion, y los tiempos de corresponden con los tiempos que tarda la onda P o S
en llegar desde cada localizacion a la estacion correspondiente.

De esta forma, si tenemos como rango de profundidades con 5 valores, los 5 primeros valores de la primera columna
se corresponderan con los tiempos de llegada de la onda desde el primer punto del mapa a la primera estacion para 
cada una de las posibles profundidades, los 5 siguientes al tiempo desde el segundo punto a cada posible 
profundidad a la primera estacion, etc...

El programa calcula los tiempos de izquierda a derecha y de arriba abajo en el mapa.

"""   
def tiempo(estaciones,puntos_mapa):
    distancias=distancia(estaciones,puntos_mapa)
    tempos=np.zeros((len(distancias[:,0])*len(puntos_mapa[2]),len(estaciones[0])))
    for k in np.arange(len(estaciones[0])):#Este bucle recorre cada estacion
        c=0
        for i in np.arange(len(distancias[:,0])): #Este buble recorre cada distancias que hemos calculado
            for j in np.arange(len(puntos_mapa[2])): #Este bucle recorre el rango de profundidades para cada distancia calculada
                #Calculamos los tiempos entre la estación k-esima y los puntos i-esimos a las profundidades j-esimas 
                a=ut.tempos(estaciones[2][k],puntos_mapa[2][j],distancias[i][k],vel,capas)
                tempos[c][k]=a[0]
                c=c+1
    return tempos

onda_P=tiempo(estaciones, puntos_mapa)
onda_S=onda_P/1.78


#Definimos la funcion mediante la que representaremos nuestra imagen
def representacion(islas_x,islas_y,est_x,est_y):
    plt.plot(islas_x,islas_y,',')     
    plt.plot(est_x,est_y,'.')
    plt.xlabel('Longitud')
    plt.ylabel('Latitud')    
    plt.grid()
    plt.show()
    
