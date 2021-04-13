#Aqui ira el programa que lea el fichero de estaciones

import matplotlib.pyplot as plt
import numpy as np
import utils_tempos as ut
import datetime as dt

#Extraemos de los ficheros las posiciones de las islas y de las esatciones sismicas
islas_long,islas_lat = np.loadtxt("contorno_islas.txt", unpack=True)
est_lat,est_long,est_alturas=np.loadtxt("fichero_estaciones.txt",comments='---',usecols=(4,5,6),unpack=True)
vel,capas=np.loadtxt("canary.txt",comments='!',usecols=(3,4),unpack=True)
nombre,fecha_i,hora_i,fecha_f,hora_f=np.genfromtxt('fichero_estaciones.txt',comments='---',usecols=(0,8,9,10,11),unpack=True,dtype=str)


rango_longitudes=np.arange(-17,-16,0.05)
rango_latitudes=np.arange(28.6,27.9,-0.05)
rango_profundidades=np.arange(0,20,5)

puntos_mapa=np.array([rango_longitudes,rango_latitudes,rango_profundidades])
est=np.array([est_long,est_lat,est_alturas])

"""
Creamos una matriz bidimensional en la que cada columna representa una profundidad.
Los datos en cada fila se ordenan de manera que se colocan para el primer punto, todos los tiempos de llegada hasta cada estacion,
despues los tiempos de llegada desde el segundo punto hasta cada esatcion y asi sucesivamente.

Si tenemos 10 estaciones, los 10 primeros valores en cada columna seran los tiempos de llegada del primer punto a cada estacion, 
y la profundidad viene indicada por la columna en la que se encuentra. Los siguientes 10 valores serían los tiempos de llegada
del segundo punto a cada una de las 10 estaciones, y asi sucesivamente.

"""


"""
Para selecionar las estaciones vamos a utilizar el modulo datetime, que nos permite comparar fechas, para 
ver cuando esta activa una estacion sismica. Primero definimos el cuando vamos a realizar el inicio de la medicion
y luego vamos a seleccionar aquellas estaciones en cuyo rango de funcionamiento entre nuestro inicio de medicion.
Para ello creamos las valiables inicio_estacion y final estacion, que son simplemente los valores temporales que 
nos dan en el fichero_estaciones pasados a formato de fecha. Lo he hecho asi (un poco raro la verdad) porque la funcion
dt.datetime.strptime no dejaba trabajar con las variables tipo str_, que son las que da por defecto python 
al extraer str de un fichero, y no lograba cambiarlas de otra manera a un str normal
"""
inicio_medicion=dt.datetime(1999,11,7,00,00,00)

def seleccion_estaciones(est,inicio_medicion):
    estaciones0=np.array([])
    estaciones1=np.array([])
    estaciones2=np.array([])
    for i in np.arange(len(est[0])):
        inicio_estacion=dt.datetime.strptime(fecha_i[i][0:10] + hora_i[i][0:7],"%Y/%m/%d%H:%M:%S")
        final_estacion=dt.datetime.strptime(fecha_f[i][0:10] + hora_i[i][0:7],"%Y/%m/%d%H:%M:%S")
        if inicio_estacion<=inicio_medicion<=final_estacion: #Comprobamos si cada i-estacion esta activa en el momento de la medicion
            estaciones0=np.append(estaciones0,est[0][i]) 
            estaciones1=np.append(estaciones1,est[0][i])
            estaciones2=np.append(estaciones2,est[0][i])  
    return estaciones0,estaciones1,estaciones2  

#He tenido que crear 3 arrays y luego juntarlos porque no me dejaba hacerlo directamente como un array de arrays 
#e ir añadiendo valores a cada uno de ellos. Tratare de hacerlo un poco mas "elegantemente mañana".


def tiempo(est,puntos_mapa):
    a,b,c=seleccion_estaciones(est,inicio_medicion)
    estaciones=np.array([a,b,c]) #Este es el array que contiene los datos geograficos de las estaciones señaladas 
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
  
#np.savetxt('Tiempos.txt',tiempo(est,puntos_mapa),fmt='%.1f',delimiter='   ')


#tiempos=tiempo(estaciones, puntos_mapa)
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
    