#Aqui ira el programa que lea el fichero de estaciones

import matplotlib.pyplot as plt
import numpy as np
import utils_tempos as ut
import   datetime as dt
import utils_hypoellipse as uh
import obspy as obs 
from scipy.interpolate import interp1d
#Extraemos de los ficheros las posiciones de las islas y de las esatciones sismicas
path = ""
islas_long,islas_lat = np.loadtxt(path +"contorno_islas.txt", unpack=True)
est_lat,est_long,est_alturas=np.loadtxt(path + "fichero_estaciones.txt",comments='---',usecols=(4,5,6),unpack=True)
vel,capas=np.loadtxt(path + "canary.txt",comments='!',usecols=(3,4),unpack=True)
nombre,fecha_i,hora_i,fecha_f,hora_f=np.genfromtxt(path + 'fichero_estaciones.txt',comments='---',usecols=(0,8,9,10,11),unpack=True,dtype=str)


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

def seleccion_estaciones(Estaciones,Momento_medicion,Estaciones_medibles,Nombre,Fecha_i,Hora_i,Fecha_f,Hora_f):
    """
    In:
        Estaciones = Array de 3 arrays en los que colocaremos la longitud, latitud y profundidad de cada estacion respectivamente
        Momento_medicion = Momento del que queremos obtener las medidas, ha de ser introducido en formato de fecha (con el modulo datetime)
        Estaciones_medibles = Array con los nombres de las estaciones que queremos tener en cuenta al realizar las medidas
        Nombre= Array con los nombres de todas las estaciones
        Fecha_i, Hora_i= Arrays con las fechas y horas del inicio de funcionamiento de cada estacion (el programa las convierte al formato datetime )
        Fecha_f, Hora_f= Arrays con las fechas y horas del final de funcionamiento de cada estacion  (el programa las convierte al formato datetime )

    Out:
        estaciones0 = array con la longitud de las estaciones seleccionadas
        estaciones1 = array con la latitud de las estaciones seleccionadas
        estaciones2 = array con la profundidad de las estaciones seleccionadas
    """
    estaciones0=np.array([])
    estaciones1=np.array([])
    estaciones2=np.array([])
    for i in np.arange(len(est[0])):
       for j in np.arange(len(estaciones_medibles)):
           if estaciones_medibles[j]==nombre[i]: #Solo vamos a usar aquellas estaciones cuyo nombre coincida con alguno de los que esta en las esatciones_medibles
               inicio_estacion=dt.datetime.strptime(fecha_i[i][0:10] + hora_i[i][0:7],"%Y/%m/%d%H:%M:%S")
               final_estacion=dt.datetime.strptime(fecha_f[i][0:10] + hora_i[i][0:7],"%Y/%m/%d%H:%M:%S")
               if inicio_estacion<=inicio_medicion<=final_estacion: #Comprobamos si cada i-estacion esta activa en el momento de la medicion
                   estaciones0=np.append(estaciones0,est[0][i]) 
                   estaciones1=np.append(estaciones1,est[1][i])
                   estaciones2=np.append(estaciones2,est[2][i])  
    return estaciones0,estaciones1,estaciones2

#He tenido que crear 3 arrays y luego juntarlos porque no me dejaba hacerlo directamente como un array de arrays 
#e ir añadiendo valores a cada uno de ellos.

def tiempo(Estaciones,Puntos_mapa,Inicio_medicion,Estaciones_medibles,Nombre,Fecha_i,Hora_i,Fecha_f,Hora_f):
    """
    In:
        Estaciones = Array de 3 arrays en los que colocaremos la longitud, latitud y profundidad[km] de cada estacion respectivamente 
        Punto_mapa = Array de 3 arrays en los que colocaremos la longitud, latitud y profundidad[km] de cada punto del mapa que queremos tener en cuenta
        Inicio_medicion = Momento del que queremos obtener las medidas, ha de ser introducido en formato de fecha (con el modulo datetime)
        Estaciones_medibles = Array con los nombres de las estaciones que queremos tener en cuenta al realizar las medidas
        Nombre= Array con los nombres de todas las estaciones
        Fecha_i, Hora_i= Arrays con las fechas y horas del inicio de funcionamiento de cada estacion (el programa las convierte al formato datetime )
        Fecha_f, Hora_f= Arrays con las fechas y horas del final de funcionamiento de cada estacion  (el programa las convierte al formato datetime )

    Out:
        tempos = Array bidimensional con los tiempos de llegada de las ondas P a las esatciones que hemos indicado que queremos usar. 
                 Cada columna representa una profundidad y en ellas se colocan los tiempos, primero desde el primer punto hasta cada 
                 una de las estaciones, despues desde el segundo punto a cada esatcion y así sucesivamente.
    """
    a,b,c=seleccion_estaciones(est,inicio_medicion,estaciones_medibles,nombre,fecha_i,hora_i,fecha_f,hora_f)
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
  
    
#Vamos a definir una funcion que concatene para cada punto dado la salida que la funcion hypoellipse_format 
def concatenador(onda_P,onda_S,estaciones_medibles,p,inicio_medicion):
    entrada_hypoellipse=""
    for i in np.arange(len(puntos_mapa[0])*len(puntos_mapa[1])): #Este bucle recorrera cada punto del mallado
        for j in np.arange(len(estaciones_medibles)): #Este bucle recorre cada estacion 	
    	    Ptime=inicio_medicion+dt.timedelta(seconds=onda_P[i+j][p]) #Añadimos para cada punto lo que tarda en llegar a cada estacion
    	    Stime=inicio_medicion+dt.timedelta(seconds=onda_S[i+j][p])
    	    a=uh.hypoellipse_format(estaciones_medibles[j],Ptime,Stime, Pw = 0, Sw = 0) #Utilizamos hypoellipse_format para que nos saque la informacion en el formato adecuado
    	    entrada_hypoellipse=entrada_hypoellipse + a + "\n" 
        entrada_hypoellipse += "\n"
    return entrada_hypoellipse


#Definimos la funcion mediante la que representaremos nuestra imagen
def representacion(islas_x,islas_y,est_x,est_y):
    plt.plot(islas_x,islas_y,',')     
    plt.plot(est_x,est_y,'.')  
    plt.xlabel('Longitud')
    plt.ylabel('Latitud')    
    plt.grid()
    plt.show()

def localizador(onda_P,onda_S,estaciones_medibles,p,inicio_medicion):
    
    """
    onda_P,onda_S: Arrays bidimensionales donde hemos almacenado los datos con todos los tiempos de llegadas de las ondas a cada estacion
    estaciones_medibles: Array de strings donde colocamos las estaciones que queremos usar para nuestras medidas.
    p: columna del array de onda P que queremos usar para las mediciones (representan cada profundidad de nuestro rango de profundidades).
    inicio_medicion: Momento del que queremos obtener las medidas, ha de ser introducido en formato de fecha (con el modulo datetime) 
    hypodata: str; String con todas las fases registradas en las estaciones que se quieran localizar
    name_file: str; Prefijo que se va a utilizar para almacenar todo
    hypoellipse_route: str; Ruta de hypoellipse
    fichero_vp_vs: str; Direccion del fichero de varaciones de la relacion Vp/Vs en funcion del tiempo, si existiese 
    hypoin_file: str; Fichero hypo*.in de configuracion de hypoellipse
    hypoctl_file: str; Fichero hypo*.ctl 
    remove: bool; Variable para guardar o no los ficheros temporales de hypoellipse, default: True
    return_values: bool; Si se quieren extraer los valores obtenidos de la localizacion configurar set True, default: True
    write_catalog_file; bool; Variable para generar catalogo propio de eventos. Actualmente solo funciona si se mete un unico evento
    magnitudes: str; String de la magnitud, solo se utiliza en el caso de localizar un unico evento. 
    
    """
      
    #Ejecutamos el concatenador para obtener la entrada del programa hypoellipse_locator y poder obtener los datos que despues debemos comparar
    hypodata=concatenador(onda_P,onda_S,estaciones_medibles,0,inicio_medicion)   
    
    #Ejecutamos la funcion hypoellipse_locator para obtener los valores que estamos buscando 
    tiempos, latitudes, longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases=\
    uh.hypoellipse_locator(hypodata, name_file = "prueba", hypoellipse_route = path + "../hypoellipse3", fichero_vp_vs = None, hypoin_file = 'hypo_hierro.in',\
    hypoctl_file = 'hypoc_hierro.ctl', remove = False,  nombre_salida = "prueba", return_values = True, write_catalog_file = False,  magnitudes = "")
    
    return tiempos, latitudes, longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases    
    

def calculo_errores(puntos_mapa,p,longitudes,latitudes,prof):
	#Ordeno las latitudes, longitudes y profundidades en un mallado 2D para restarlo mas facilmente con las obtenidas por hypoellipse
	lat_teoricas=np.zeros((len(puntos_mapa[1]),len(puntos_mapa[0])))
	long_teoricas=np.zeros((len(puntos_mapa[1]),len(puntos_mapa[0])))
	prof_teoricas=np.zeros((len(puntos_mapa[1]),len(puntos_mapa[0])))+p
	for i in np.arange(len(puntos_mapa[1])):
		for j in np.arange(len(puntos_mapa[0])):
			long_teoricas[i][j]=puntos_mapa[0][j]
			lat_teoricas[i][j]=puntos_mapa[1][i]
	long_hypo=np.array(longitudes).reshape(len(puntos_mapa[1]),len(puntos_mapa[0]))			
	lat_hypo=np.array(latitudes).reshape(len(puntos_mapa[1]),len(puntos_mapa[0]))
	prof_hypo=np.array(prof).reshape(len(puntos_mapa[1]),len(puntos_mapa[0]))
	errores_lat=abs(lat_teoricas-lat_hypo)*111
	errores_long=abs(long_teoricas-long_hypo)*111
	errores_prof=abs(prof_teoricas-prof_hypo)
	return errores_lat,errores_long,errores_prof

#Definimos algunas de nuestras variables
inicio_medicion=dt.datetime(2021,1,1,00,00,00)
estaciones_medibles=np.array(["CCAN","CCHO","CGRA","CLAJ","CPVI", "CTFS","MACI"])   
onda_P=tiempo(est,puntos_mapa,inicio_medicion,estaciones_medibles,nombre,fecha_i,hora_i,fecha_f,hora_f)
onda_S=onda_P*1.78
 
#Ejecutamos la funcion localizador para obtener los ficheros que necesitamos    
tiempos, latitudes,longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases=localizador(onda_P,onda_S,estaciones_medibles,0,inicio_medicion)

#Transformamos las latitudes, longitudes y profundidades a floats, ya que el programa las saca como strings 
longitudes=longitudes.astype(float)
latitudes=latitudes.astype(float)
prof=prof.astype(float)

#ejecutamos la funcion mediante la que calculamos los errores
errores_lat,errores_long,errores_prof=calculo_errores(puntos_mapa,0,longitudes,latitudes,prof)			

"""	
np.savetxt('error_long', errores_long,fmt='%.2f')
np.savetxt('error_lat', errores_lat,fmt='%.2f')
np.savetxt('error_prof', errores_prof,fmt='%.2f')
"""    
    
    
    
    
    
    
