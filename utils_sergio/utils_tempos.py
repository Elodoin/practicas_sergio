import numpy as np 
from scipy.optimize import brentq


# Calculo de los desfases temporales 
def tempos(H, h, D, v=np.array([4.2,6.3,7.5,8,8,8]),  z=np.array([0,4,12,18,20,32]) ):
    #Redondeamos para evitar problemas con brentq() despues.
    h = np.around(h, 10)
    H = np.around(H, 10)
    D = np.around(D, 10)
    #Alturas y velocidades de capas en Canarias (en km y km/s). Trabajaremos con la velocidad de la onda p por defecto, pudiendo hallar el tiempo de la s multiplicando la salida por 1.78.


    velocidad = 0
    if h > 0:
        velocidad = v[ np.max( np.arange( len(z) )[(z < h) | (z == h)] ) ]
    if h < 0:
        velocidad = v[0]
    #En el array tempos guardaremos los tiempos de llegada de cada rayo y en el array capa la capa en que se produce la refraccion critica, asignandole el valor de 0 al rayo directo.
    tempos = []
    pes = []
    capa  = []
    
    #Distinguimos ahora los casos en que el terremoto se produce a una altura superior a la del sismometro y viceversa. Comenzamos con el caso en que el terremoto se produce por debajo del sismometro:
    if h + H >= 0:
        #Situamos el origen de las capas en el sismometro, de forma que variamos la profundidad de las capas y llamamos h a la profundidad del terremoto para este nuevo sistema de referencia.
        z = np.concatenate(([0], z[1:] + H))

	#La siguiente linea de codigo nos solventa el problema que surge cuando el sismometro se produce a una profundidad > a z[1].
        z1 = z[z>=0]
        #n0 representa el indice original de la primera capa.
        n0 = len(z)-len(z1)
        v = v[n0:]
        #Redefinimos h para el nuevo sistema de referencia.
        h = h + H

        #Separamos las capas entre aquellas por encima y por debajo del evento.
        direc = z1[z1<h]
        refrac = z1[z1>h]

        #Para estudiar el rayo directo, a las que estan por encima le anadimos la profundidad del evento y calculamos las distancias entre capas.
        direc1 = np.concatenate((direc, [h]))
        zp = direc1[1:]-direc1[:-1]

        if h > 0 :
            #Introducimos una funcion que nos permite hallar el parametro p (equivalente a hallar el angulo de incidencia) para el rayo directo.
            def func(p):
                u = D - np.sum(zp/np.sqrt((v[:len(zp)]*p)**(-2)-1))
                return u

            #Buscamos el cero de esta funcion mediante el metodo numerico brentq. (Ojo!, el metodo converge bien pero en ciertos casos advierte de una division por cero que puede ralentizar el programa)
            p = brentq(func, 0, (1./(np.max(v[:len(zp)]))))

            #Hallamos el tiempo que tarda en llegar el rayo directo.
            t = np.sum(zp/(v[:len(zp)]*np.sqrt(1-(p*v[:len(zp)])**2)))

        #Consideramos aqui el caso del rayo horizontal.
        else:
            direc = np.array([0])
            t = D/v[0]
            p = 0.

        #Guardamos los resultados.
        tempos.append(t)
        capa.append(0)
        pes.append(p)

        #Vamos a estudiar ahora el tiempo empleado por los diferentes rayos refractados.
        for n in range(1, len(refrac)+1):
            #En direc2 guardamos las capas que atraviesa el rayo cuando sube desde la refraccion critica hasta el sismometro (trayectoria de subida); y en refrac2 a las capas que atraviesa al bajar desde su hipocentro hasta la refraccion critica (trayectoria de bajada). Ademas, a la trayectoria de bajada le anadimos la profundidad del terremoto. Hallamos a continuacion las distancias entre capas para ambas trayectorias.
            direc2 = np.concatenate((direc, refrac[:n]))
            refrac2 = refrac[:n]
            zp = direc2[1:]-direc2[:-1]
            refrac2 = np.concatenate(([h], refrac2))
            zr = refrac2[1:]-refrac2[:-1]

            #Se puede comprobar que para estos rayos que sufren refraccion critica el valor del parametro p viene dado por p = 1/vn, siendo vn la velocidad de propagacion en la capa en que se produce la refraccion critica.
            p = 1/v[len(direc2)-1]

            #Introducimos una funcion que nos permite hallar la distancia dn que recorre el rayo por la capa de refraccion critica.
            def fanc(dn):
                u = D - np.sum(zp/np.sqrt((v[:len(zp)]*p)**(-2)-1)) - np.sum(zr/np.sqrt((v[len(direc)-1: n+len(direc)-1]*p)**(-2)-1)) - dn
                return u

            #Descartamos aquellas refracciones que no pueden darse pues se supera el angulo critico: fanc(0)<= 0:
            if fanc(0) > 0:
                #Buscamos numericamente el valor de dn para cada rayo y hallamos el tiempo empleado, donde tenemos en cuanto las tres trayectorias distinguibles seguidas por el rayo: bajada, critica y subida, guardandolo en la lista tempo asi como la capa critica que alcanzo.
                dn = brentq(fanc, 0, D)
                t = np.sum(zp/(v[:len(zp)]*np.sqrt(1-(p*v[:len(zp)])**2))) + np.sum(zr/(v[len(direc)-1: n+len(direc)-1]*np.sqrt(1-(p*v[len(direc)-1: n+len(direc)-1])**2))) + dn*p
                tempos.append(t)
                capa.append(len(direc) + n -1)
                pes.append(p)
                #El numero guardado corresponde con la capa critica: 2 -> z2, que es la segunda capa (cuidado con los indices!)

            else:
                break

    #Vamos a estudiar ahora el caso en que el seismo se produce por encima del sismometro.
    else:
        #Cambiamos el origen de coordenadas de forma que ahora el terrremoto se produzca a una altitud 0, y el sismometro se encuentre a una profundidad h-H.
        z = np.concatenate(([0], z[1:] - h))
        z1 = z[ z >= 0]
        #n0 representa el indice original de la primera capa.
        n0 = len(z) - len(z1)
        v = v[n0:]
        H = abs(h + H)

        #Separamos las capas entre aquellas por encima y por debajo del evento.
        direc = z1[z1 < H]
        refrac = z1[z1 > H]

        #Para estudiar el rayo directo, a las que estan por encima le anadimos la profundidad del evento y calculamos las distancias entre capas.
        direc1 = np.concatenate((direc, [H]))
        zp = direc1[1:] - direc1[:-1]

        #Introducimos una funcion que nos permite hallar el parametro p (equivalente a hallar el angulo de incidencia) para el rayo directo.
        def func(p):
            u = D - np.sum(zp / np.sqrt((v[:len(zp)]*p)**(-2)-1))
            return u

        #Buscamos el cero de esta funcion mediante el metodo numerico brentq. (Ojo!, el metodo converge bien pero en ciertos casos advierte de una division por cero que puede ralentizar el programa)
        p = brentq(func, 0, (1 / (max(v[:len(zp)]))))

        #Hallamos el tiempo que tarda en llegar el rayo directo.
        t = np.sum(zp / (v[:len(zp)]*np.sqrt(1 - (p*v[:len(zp)])**2)))

        tempos.append(t)
        capa.append(0)
        pes.append(p)

        #Para las refractadas seguimos un procedimiento similar al caso anterior.
        for n in np.arange(1, len(refrac)+1):
            #Al contrario que antes, ahora la direc2 constituye el camino de bajada y la refrac2 el de subida, de resto todo es analogo.
            direc2 = np.concatenate((direc, refrac[:n]))
            refrac2 = refrac[:n]
            zp = direc2[1:] - direc2[:-1]
            refrac2 = np.concatenate(([H], refrac2))
            zr = refrac2[1:] - refrac2[:-1]

            #Al igual que antes, el valor de p viene dado por la capa de la refraccion critica.
            p = 1 / v[len(direc2) - 1]

            #Introducimos una funcion que nos permite hallar la distancia dn que recorre el rayo por la capa de refraccion critica.
            def fanc(dn):
                u = D - np.sum(zp/np.sqrt((v[:len(zp)]*p)**(-2)-1)) - np.sum(zr/np.sqrt((v[len(direc)-1: n+len(direc)-1]*p)**(-2)-1)) - dn
                return u

            #Descartamos aquellas refracciones que no pueden darse pues se supera el angulo critico: fanc(0)<= 0:
            if fanc(0) > 0:
                #Buscamos numericamente el valor de dn para cada rayo y hallamos el tiempo empleado, donde tenemos en cuanto las tres trayectorias distinguibles seguidas por el rayo: bajada, critica y subida, guardandolo en la lista tempo asi como la capa critica que alcanzo.
                dn = brentq(fanc, 0, D)
                t = np.sum(zp/(v[:len(zp)]*np.sqrt(1-(p*v[:len(zp)])**2))) + np.sum(zr/(v[len(direc)-1: n+len(direc)-1]*np.sqrt(1-(p*v[len(direc)-1: n+len(direc)-1])**2))) + dn*p
                tempos.append(t)
                capa.append(len(direc) + n - 1)
                pes.append(p)
                #El numero guardado corresponde con la capa critica: 2 -> z2, que es la segunda capa (cuidado con los indices!)

            else:
                break

    return np.array([min(tempos)])[0], velocidad, pes[np.argmin(tempos)]


# Calculo de la distancia entre dos puntos en la superficie de la tierra 
def distance(lat1, lon1, lat2, lon2):
    # Earth Radius in KM
    R = 6371
    dLat = np.radians(lat2 - lat1)
    dLon = np.radians(lon2 - lon1)
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    a = np.sin(dLat / 2) * np.sin(dLat / 2) +  np.sin(dLon / 2) * np.sin(dLon / 2) * np.cos(lat1) * np.cos(lat2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = R * c
    return d


