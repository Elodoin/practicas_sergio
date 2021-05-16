#Aqui ira el programa que lea el fichero de estaciones

import matplotlib.pyplot as plt
import numpy as np
import utils_tempos as ut
import   datetime as dt
import utils_hypoellipse as uh
import obspy as obs 
from scipy.interpolate import interp1d
import random as ra
#Extraemos de los ficheros las posiciones de las islas y de las esatciones sismicas
path = ""
islas_long,islas_lat = np.loadtxt(path +"contorno_islas.txt", unpack=True)
est_lat,est_long,est_alturas=np.loadtxt(path + "stations.txt",comments=('---'),usecols=(4,5,6),unpack=True)
vel,capas=np.loadtxt(path + "canary.txt",comments='!',usecols=(3,4),unpack=True)
nombre,fecha_i,hora_i,fecha_f,hora_f=np.genfromtxt(path + 'stations.txt',comments='---',usecols=(0,8,9,10,11),unpack=True,dtype=str)

rango_longitudes=np.arange(-17.1,-15.9,0.05)
rango_latitudes=np.arange(28.8,27.9,-0.05)
rango_profundidades=np.arange(0,20,15)

puntos_mapa=np.array([rango_longitudes,rango_latitudes,rango_profundidades])
est=np.array([est_long,est_lat,est_alturas])

#Definimos algunas de nuestras variables
inicio_medicion=dt.datetime(2021,1,1,00,00,00)
estaciones_entrada=np.sort(list(nombre[0:30]))


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

def seleccion_estaciones(Estaciones,Momento_medicion,Estaciones_entrada,Nombre,Fecha_i,Hora_i,Fecha_f,Hora_f):
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
    estaciones_medibles=[]
    #Creamos un array con los nombres de las estaciones que queremos medir sin que se repita ningun nombre
    for item in Nombre:
    	if item not in estaciones_medibles: 
    		if item in Estaciones_entrada:
        		estaciones_medibles=np.append(estaciones_medibles,item) #Estaciones cuyo nombre coincide con las que damos por la entrada

    i=0
    for n in Nombre:
    	
    	if n in estaciones_medibles:
    	 	inicio_estacion=dt.datetime.strptime(fecha_i[i][0:10] + hora_i[i][0:7],"%Y/%m/%d%H:%M:%S")
    	 	final_estacion=dt.datetime.strptime(fecha_f[i][0:10] + hora_i[i][0:7],"%Y/%m/%d%H:%M:%S")
    	 	if inicio_estacion<=inicio_medicion<=final_estacion: #Comprobamos si cada i-estacion esta activa en el momento de la medicion
    	 		estaciones0=np.append(estaciones0,est[0][i])
    	 		estaciones1=np.append(estaciones1,est[1][i])
    	 		estaciones2=np.append(estaciones2,est[2][i])
    	i=i+1
    return estaciones0,estaciones1,estaciones2,estaciones_medibles

#He tenido que crear 3 arrays y luego juntarlos porque no me dejaba hacerlo directamente como un array de arrays 
#e ir añadiendo valores a cada uno de ellos.

def extractor_residuos(nombre_fichero_out,estaciones_medibles):
	a,b,c,estaciones_medibles=seleccion_estaciones(est,inicio_medicion,estaciones_entrada,nombre,fecha_i,hora_i,fecha_f,hora_f)
		
	fichero=open("../hypoellipse3/" + "%s" % nombre_fichero_out, 'r').readlines()	
	busqueda='stn c pha remk p p-sec s-sec resid  std-er   dist  azm ain    tc c vthk  ttob-ttcal-dlay-edly=resid rmk stn pha sources'
	l=0
	indices=[]
	residuos_S=[]
	residuos_P=[]
	for linea in fichero:
		if busqueda in linea:
			indices=np.append(indices,l)
		l=l+1
	indices=indices.astype(int)	
	for i in indices:
		datos_estaciones=fichero[i+1:i+1+2*len(estaciones_medibles)]
		
		for j in range(len(estaciones_medibles)):
			n=0
			for fila in datos_estaciones:
			
				if estaciones_medibles[j].swapcase() in fila and n%2==0:
				
					residuos_P=np.append(residuos_P,float(fila[32:36]))
					 
				if estaciones_medibles[j].swapcase() in fila and n%2!=0:
					
					residuos_S=np.append(residuos_S,float(fila[32:36]))
				n+=1
	return residuos_P,residuos_S
	







def tiempo(Estaciones,Puntos_mapa,Inicio_medicion,estaciones_entrada,Nombre,Fecha_i,Hora_i,Fecha_f,Hora_f):
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
                 una de las estaciones, despues desde el segundo punto a cada estacion y así sucesivamente.
    """
    a,b,c,estaciones_medibles=seleccion_estaciones(est,inicio_medicion,estaciones_entrada,nombre,fecha_i,hora_i,fecha_f,hora_f)
    estaciones=np.array([a,b,c]) #Este es el array que contiene los datos geograficos de las estaciones señaladas 
    d=np.zeros((len(estaciones[0]),len(puntos_mapa[0])*len(puntos_mapa[1]))) #Matriz donde cada columna es un punto del mapa y cada fila una estacion
    tempos=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1])*len(estaciones_medibles),len(puntos_mapa[2])))
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

    return tempos,estaciones,d,estaciones_medibles
  
    
#Vamos a definir una funcion que concatene para cada punto dado la salida que la funcion hypoellipse_format 
def concatenador(onda_P,onda_S,p,inicio_medicion,estaciones_medibles):
    entrada_hypoellipse=""
    b=0
    for i in np.arange(len(puntos_mapa[0])*len(puntos_mapa[1])): #Este bucle recorrera cada punto del mallado
        for j in np.arange(len(estaciones_medibles)): #Este bucle recorre cada estacion 	
    	    Ptime=inicio_medicion+dt.timedelta(seconds=onda_P[b][p]) #Añadimos para cada punto lo que tarda en llegar a cada estacion
    	    Stime=inicio_medicion+dt.timedelta(seconds=onda_S[b][p])
    	    a=uh.hypoellipse_format(estaciones_medibles[j],Ptime,Stime, Pw = 0, Sw = 0) #Utilizamos hypoellipse_format para que nos saque la informacion en el formato adecuado
    	    entrada_hypoellipse=entrada_hypoellipse + a + "\n"

    	    b=b+1
    	    
        entrada_hypoellipse += "\n"
    return entrada_hypoellipse




def localizador(onda_P,onda_S,p,inicio_medicion,estaciones_medibles):
    
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
    hypodata=concatenador(onda_P,onda_S,p,inicio_medicion,estaciones_medibles)   
    #Ejecutamos la funcion hypoellipse_locator para obtener los valores que estamos buscando 
    tiempos, latitudes, longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases=\
    uh.hypoellipse_locator(hypodata, name_file = "prueba", hypoellipse_route = path + "../hypoellipse3", fichero_vp_vs = None, hypoin_file = 'hypo_hierro.in',\
    hypoctl_file = 'hypoc_hierro.ctl', remove = False,  nombre_salida = "prueba", return_values = True, write_catalog_file = False,  magnitudes = "")
    
    r_P,r_S=extractor_residuos("prueba_tmp.out",estaciones_medibles)
    
    
    return tiempos, latitudes, longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases,r_P,r_S    
    

def aplicacion_residuos(est,puntos_mapa,inicio_medicion,estaciones_entrada,nombre,fecha_i,hora_i,fecha_f,hora_f):
	
	#Obtenemos los valores de llegada de las ondas P y S sin aplicarles residuos
    	onda_P,estaciones,d,estaciones_medibles=tiempo(est,puntos_mapa,inicio_medicion,estaciones_entrada,nombre,fecha_i,hora_i,fecha_f,hora_f)
    	onda_S=onda_P*1.78
    	
    	#Generamos los arrays donde almacenaremos los residuos
    	residuos_P=np.zeros((len(onda_P),len(puntos_mapa[2])))
    	residuos_S=np.zeros((len(onda_P),len(puntos_mapa[2])))

    	#Definimos el numero de veces que queremos calcular las coordenadas con los residuos aplicados aleatoriamente
    	n=5
    	
    	#Creamos los arrays en los que meteremos los distintos valores de las latitudes y longitudes tras sumarles aleatoriamente los residuos a las ondas P y S
    	
    	lats=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1]),n*len(puntos_mapa[2])))
    	longs=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1]),n*len(puntos_mapa[2])))
    	profs=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1]),n*len(puntos_mapa[2])))
    	#Generamos los residuos
    	for p in range(len(puntos_mapa[2])):
    		tiempos, latitudes,longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap,numero_fases,r_P,r_S=\
    		localizador(onda_P,onda_S,p,inicio_medicion,estaciones_medibles)
    		#Colocamos los residuos en los arrays que habiamos definido
    		residuos_P[:,p]=r_P
    		residuos_S[:,p]=r_S
    	a=0
    	#Generamos la matriz 2D en la que en las primeras p-columnas estaran los datos de la 1º aplicacion de los residuos para cada profundidad p
    	#En las siguientes p-columnas estaran los datos de la 2º aplicacion de los residuos de manera aleatoria para cada profundidad p, y asi sucesivamente
    	for i in range(n):
    		
    		for j in range(len(residuos_P[:,0])):
    			for k in range(len(residuos_P[0,:])):
    				residuos_P[j][k]=residuos_P[j][k]*ra.choice((1,-1))
    				residuos_S[j][k]=residuos_S[j][k]*ra.choice((1,-1))

    		for p in range(len(puntos_mapa[2])):
    			o_P=onda_P+residuos_P
    			o_S=onda_S+residuos_P
    			tiempos, latitudes,longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap,numero_fases,r_P,r_S=\
    			localizador(o_P,o_S,p,inicio_medicion,estaciones_medibles)
    			lats[:,a]=latitudes
    			longs[:,a]=longitudes
    			profs[:,a]=prof
    			a+=1
    	media_lat=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1]),len(puntos_mapa[2])))
    	media_long=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1]),len(puntos_mapa[2])))
    	media_prof=np.zeros((len(puntos_mapa[0])*len(puntos_mapa[1]),len(puntos_mapa[2])))	
    	#Asi que tenemos que hacer una media entre los valores de las columnas p,p+len(puntos_mapa[2]),p+2*len(puntos_mapa[2]),....,p+n*len(puntos_mapa[2])
    	for i in range(len(lats[:,0])):
    		for p in range(len(puntos_mapa[2])):
    			la=[]
    			lo=[]
    			pr=[]
    			for l in range(n):
    				la=np.append(la,lats[i,p+l*len(puntos_mapa[2])])
    				lo=np.append(lo,longs[i,p+l*len(puntos_mapa[2])])
    				pr=np.append(pr,profs[i,p+l*len(puntos_mapa[2])])
    			media_lat[i,p]=np.mean(la)
    			media_long[i,p]=np.mean(lo)
    			media_prof[i,p]=np.mean(pr)
    			
    	#Transformamos las latitudes, longitudes y profundidades a floats, ya que el programa las saca como strings y las coloco en un mallado 2D
    	media_lat=media_lat.astype(float)
    	media_long=media_long.astype(float)
    	media_prof=media_prof.astype(float)    
    			
    	return media_lat,media_long,media_prof,estaciones


def calculo_errores(p,est,puntos_mapa,inicio_medicion,estaciones_entrada,nombre,fecha_i,hora_i,fecha_f,hora_f):
	
	media_lat,media_long,media_prof,estaciones=aplicacion_residuos(est,puntos_mapa,inicio_medicion,estaciones_entrada,nombre,fecha_i,hora_i,fecha_f,hora_f)
	long_media_hypo=np.array(media_long[:,p]).reshape(len(puntos_mapa[0]),len(puntos_mapa[1])).T
	lat_media_hypo=np.array(media_lat[:,p]).reshape(len(puntos_mapa[0]),len(puntos_mapa[1])).T
	prof_media_hypo=np.array(media_prof[:,p]).reshape(len(puntos_mapa[0]),len(puntos_mapa[1])).T
	
	#Ordeno las latitudes, longitudes y profundidades teoricas en un mallado 2D para restarlo mas facilmente con las obtenidas por hypoellipse
	lat_teoricas=np.zeros((len(puntos_mapa[1]),len(puntos_mapa[0])))
	long_teoricas=np.zeros((len(puntos_mapa[1]),len(puntos_mapa[0])))
	prof_teoricas=np.zeros((len(puntos_mapa[1]),len(puntos_mapa[0])))+puntos_mapa[2][p]
	
	for i in np.arange(len(puntos_mapa[1])):
		for j in np.arange(len(puntos_mapa[0])):
			long_teoricas[i][j]=puntos_mapa[0][j]
			lat_teoricas[i][j]=puntos_mapa[1][i]
	errores_lat=abs(lat_teoricas-lat_media_hypo)*111
	errores_long=abs((long_teoricas-long_media_hypo)*np.cos(lat_teoricas*(np.pi/180.))*111)
	errores_prof=abs(prof_teoricas-prof_media_hypo)

	return errores_lat,errores_long,errores_prof,estaciones

#Definimos la funcion mediante la que representaremos nuestra imagen
def representacion(islas_long,islas_lat,est,puntos_mapa,estaciones_entrada,inicio_medicion,nombre,fecha_i,hora_i,fecha_f,hora_f):
	    
    			
    	nombre_errores=['longitud','latitud','profundidad']	

    	fig, axs = plt.subplots(len(puntos_mapa[2]),3)
	
    	for p in range(len(puntos_mapa[2])):	
    		errores_lat,errores_long,errores_prof,estaciones=calculo_errores(p,est,puntos_mapa,inicio_medicion,estaciones_entrada,nombre,fecha_i,hora_i,fecha_f,hora_f)		
    		errores=np.array([errores_long,errores_lat,errores_prof])
    		for col in range(3):
    	    		if len(puntos_mapa[2])==1:
    	    			ax=axs[col]
    	    		else:
    	    			ax= axs[p,col]
    	    		a=ax.imshow(errores[col],extent=[min(puntos_mapa[0]),max(puntos_mapa[0]),min(puntos_mapa[1]),max(puntos_mapa[1])],alpha=0.5,cmap='viridis')
    	    		ax.plot(islas_long,islas_lat,color='k')
    	    		ax.plot(estaciones[0],estaciones[1],'.',color='g')
    	    		ax.set_xlabel('Longitud')
    	    		ax.set_ylabel('Latitud')
    	    		ax.set_title('Error en %s a profundidad %i [km]' %(nombre_errores[col],puntos_mapa[2][p]))    
    	    		ax.set_xlim(min(puntos_mapa[0]),max(puntos_mapa[0]))
    	    		ax.set_ylim(min(puntos_mapa[1]),max(puntos_mapa[1])) 
    	    		fig.colorbar(a,ax=ax,shrink=0.4)
    	    		fig.tight_layout()
    	    		
    	plt.show()

#Ejecutamos el programa de representacion de los mapas de error

representacion(islas_long,islas_lat,est,puntos_mapa,estaciones_entrada,inicio_medicion,nombre,fecha_i,hora_i,fecha_f,hora_f)


    
