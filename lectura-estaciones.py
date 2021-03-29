#Aqui ira el programa que lea el fichero de estaciones

import matplotlib.pyplot as plt
import numpy as np
import pandas as pa

#Extraemos de los ficheros las posiciones de las islas y de las esatciones sismicas
islas_x,islas_y = np.loadtxt("utils_sergio//contorno_islas.txt", unpack=True)
est_y,est_x=np.loadtxt("utils_sergio//fichero_estaciones.txt",comments='---',usecols=(4,5),unpack=True)

#Definimos la funcion mediante la que representaremos nuestra imagen
def representacion(islas_x,islas_y,est_x,est_y):
    plt.plot(islas_x,islas_y,',')     
    plt.plot(est_x,est_y,'.')
    plt.xlabel('Longitud')
    plt.ylabel('Latitud')    
    plt.grid()
    plt.show()
    
representacion(islas_x,islas_y,est_x,est_y)    