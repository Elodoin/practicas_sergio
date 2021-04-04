# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 19:36:04 2021

@author: Sergio Catalán
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pa
import utils_tempos as ut
import random
#import cv2
 
# creamos una matriz de 100x100 con todos los colores en blanco
img = np.zeros((100,100,3),np.uint8)
 
# recorremos cada uno de los elementos
for x in range(100):
    for y in range(100):
 
        # Cambiamos el color de cada uno de los pixeles de forma aleatoria
        img[x,y] = [random.randint(0,256),random.randint(0,256),random.randint(0,256)]

plt.plot(img)
# guardamos la imagen png
#cv2.imwrite("MyImage.png",img)




"""

#Extraemos de los ficheros las posiciones de las islas y de las estaciones sismicas
islas_x,islas_y = np.loadtxt("contorno_islas.txt", unpack=True)
est_y,est_x=np.loadtxt("fichero_estaciones.txt",comments='---',usecols=(4,5),unpack=True)
#Tambien las velocidades y profundidades de cada capa
vel,capas=np.loadtxt("canary.prm",comments='!',usecols=(3,4),unpack=True)

#ut.distance()


#Definimos la funcion mediante la que representaremos nuestra imagen
def representacion(islas_x,islas_y,est_x,est_y):
    plt.plot(islas_x,islas_y,',')     
    plt.plot(est_x,est_y,'.')
    plt.xlabel('Longitud')
    plt.ylabel('Latitud')    
    plt.grid()
    plt.show()
    
representacion(islas_x,islas_y,est_x,est_y)    
"""