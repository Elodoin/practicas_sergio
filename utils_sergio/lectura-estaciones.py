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
rango_profundidades=np.arange(0,20,5)

puntos_mapa=np.array([rango_longitudes,rango_latitudes,rango_profundidades])
estaciones=np.array([est_long,est_lat,est_alturas])

"""
Creamos una matriz bidimensional en la que cada columna representa una profundidad.
Los datos en cada fila se ordenan de manera que se colocan para el primer punto, todos los tiempos de llegada hasta cada estacion,
despues los tiempos de llegada desde el segundo punto hasta cada esatcion y asi sucesivamente.

Si tenemos 10 estaciones, los 10 primeros valores en cada columna seran los tiempos de llegada del primer punto a cada estacion, 
y la profundidad viene indicada por la columna en la que se encuentra. Los siguientes 10 valores serían los tiempos de llegada
del segundo punto a cada una de las 10 estaciones, y asi sucesivamente.
"""

def tiempo(estaciones,puntos_mapa):
    d=np.zeros((len(estaciones[0]),len(puntos_mapa[0])*len(puntos_mapa[1]))) #Matriz donde cada columna es un punto del mapa y cada fila una estacion
    tempos=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1])*len(estaciones[0]),len(puntos_mapa[2])))
    c=0 #la c representa cada punto del mapa
    p=0 #la p representa cada fila de datos
    for i in np.arange(len(puntos_mapa[0])): #Este buble recorre cada longitud en el mapa (eje x)
        for j in np.arange(len(puntos_mapa[1])): #Este bucle recorre cada latitud en el mapa (eje y)
            for r in np.arange(len(estaciones[0])): #Este bucle recorre cada estacion
                d[r][c]=ut.distance(puntos_mapa[1][j],puntos_mapa[0][i],estaciones[1][r],estaciones[0][r]) #Calculamos las distancias
                for k in np.arange(len(puntos_mapa[2])): #Este bucle recorre cada profundidad
                    t=ut.tempos(estaciones[2][r]/1000,puntos_mapa[2][k],d[r][c],vel,capas) 
                    tempos[p][k]=t[0]
                p+=1 
            c+=1    
    return tempos
  


tiempos=tiempo(estaciones, puntos_mapa)
#onda_P=tiempo(estaciones, puntos_mapa)
#onda_S=onda_P/1.78
#np.savetxt('OndaP.txt',onda_P,fmt='%.1f',delimiter='   ')

#Definimos la funcion mediante la que representaremos nuestra imagen
def representacion(islas_x,islas_y,est_x,est_y):
    plt.plot(islas_x,islas_y,',')     
    plt.plot(est_x,est_y,'.')
    plt.xlabel('Longitud')
    plt.ylabel('Latitud')    
    plt.grid()
    plt.show()
    