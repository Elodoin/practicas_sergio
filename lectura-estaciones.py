#Aqui ira el programa que lea el fichero de estaciones


import pylab as py
islas_x,islas_y = py.loadtxt("utils_sergio\\contorno_islas.txt", unpack=True)
est_y,est_x=py.loadtxt("utils_sergio\\fichero_estaciones.txt",comments='---',usecols=(4,5),unpack=True)

py.ion()
py.plot(islas_x,islas_y,',')     
py.plot(est_x,est_y,'.')
py.xlabel('Longitud')
py.ylabel('Latitud')    
py.grid()
py.show()