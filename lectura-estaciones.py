#Aqui ira el programa que lea el fichero de estaciones


import pylab as py
islas_x,islas_y = py.loadtxt("C:\\Users\\Sergio Catalán\\Documents\\GitHub\\practicas_sergio\\utils_sergio\\contorno_islas.txt", unpack=True)
est_y,est_x=py.loadtxt("C:\\Users\\Sergio Catalán\\Documents\\GitHub\\practicas_sergio\\utils_sergio\\fichero_estaciones.txt",comments='---',usecols=(3,4),unpack=True)


py.ion()
py.plot(islas_x,islas_y,',')     
py.plot(est_x,est_y,'.','r')    
py.xlim(min(islas_x-0.3),max(islas_x+0.3))
py.ylim(min(islas_y-0.3),max(islas_y+0.3))
py.show()