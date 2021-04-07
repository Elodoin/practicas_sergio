import os
import numpy as np
import obspy as obs 
from scipy.interpolate import interp1d


# Para transformar los picados a valores de hypoellipse_format
def hypoellipse_data(Values, Stations):

    events = []
    keys0 = np.array([*Values.keys()])
    Stations = np.array([" "*(4-len(sta.replace(" ", ""))) + sta for sta in Stations])
    keys = keys0[np.argsort(Stations[keys0])]

    for key in keys:
        StatDict = Values[key]
        Pvalue, Svalue = 0,0
        Pw, Sw = StatDict[3]

        if StatDict[2][0] == 0:
            Pw=5

        if StatDict[2][1] == 0:
            Sw=5

        Pvalue = StatDict[2][0]
        Svalue = StatDict[2][1]
        events.append(hypoellipse_format(Stations[key].replace(" ", "").replace("E020", "E02"), obs.UTCDateTime(Pvalue),obs.UTCDateTime( Svalue), Pw, Sw) + "\n")
 
    events.sort()

    events.append("\n")
    u = ""
    for rstr in events:
        u=u+rstr
    return u


def hypoellipse_data_ign(fecha, event, stations_to_use = str):
    
    phrase = "" 
    
    stations = np.sort(list(set(event.T[0])))
    
    stations_to_compare = np.array([sta.replace(" ", "") for sta in stations])
    
    if type(stations_to_use) != str: 
        stations = stations[np.in1d(stations_to_compare, stations_to_use)]
    
    stations = stations[np.argsort([" "*(4-len(sta.replace(" ", "")))+sta.replace(" ", "") for sta in stations])]
    
    for sta in stations:
        Ptime = obs.UTCDateTime(0)
        Stime = obs.UTCDateTime(0)
        Pw = 5
        Sw = 5
        station_values = event[event.T[0] == sta]
        p_masked = station_values[station_values[:, 3] == "P"]
        s_masked = station_values[(station_values[:, 3] == "S") | (station_values[:, 3] == "L")]

        if len(p_masked) != 0:
            Pw = 0
            Ptime = obs.UTCDateTime(fecha + " " + p_masked[0][4])

        if len(s_masked) != 0:
            Sw = 0
            Stime = obs.UTCDateTime(fecha + " " + s_masked[0][4])
       
       
        phrase += hypoellipse_format(sta.replace(" ", "").replace("E020", "E02"), Ptime, Stime, Pw, Sw) + "\n"

    return phrase + "\n"


# Escritura de fases en version hypoellipse
def hypoellipse_format(Station, Ptime, Stime, Pw = 0, Sw = 0):
    """
    Station: str; Nombre de la estacion
    Ptime: datetime; Tiempo de la onda P
    Stime: datetime; Tiempo de la onda S 
    Pw: int; Peso de la onda P 
    Sw: int; Peso de la onda s
    """
    # Como arreglar el caso en el que hay un picado
    flag = 2

    if (Pw == 5) & (Sw < 5):
        Pw = " "
        Ptime = Stime
        flag = 0

    elif (Pw < 5) & (Sw == 5):
        Sw = " "
        Stime = Ptime
        flag = 1

    Pval = Ptime.datetime
    Sval = Stime.datetime

	# Ahora introducimos defectos de forma del formato de hypoellipse
    year = str(Pval.year)[2:]
    day = str(Pval.day)
    month = str(Pval.month)
    hour = str(Pval.hour)
    minu = str(Pval.minute)
    sec = str(Pval.second)

    if Pval.day < 10:
    	day  = "0" + str(Pval.day)
    if Pval.month < 10:
        month = "0" + str(Pval.month)
    if Pval.hour < 10:
        hour = "0" + str(Pval.hour)
    if Pval.minute < 10:
        minu = "0" + str(Pval.minute)
    pstr = year + month + day + hour + minu + " "
    #Asumimos en el truncamiento siguiente un error de 0.005 s, lo cual es necesario pues hypoellipse en su entrada solo acepta hasta el centisegundo.
    if flag != 0:
        isec = str(int(Pval.microsecond/10e3))
        if int(isec)<10:
            isec = "0" + isec
        pr = sec + isec
        if len(pr)<4:
            pr = " " +  pr

    if flag == 0:
        pr = "    "

    #Asumimos en el truncamiento siguiente un error de 0.005 s, lo cual es necesario pues hypoellipse en su entrada solo acepta hasta el centisegundo.
    if flag !=1:
        
        ssec = str(Sval.second)
        sisec = str(int(Sval.microsecond/10e3))
        if int(sisec)<10:
            sisec = "0" + sisec
        sr = ssec + sisec
        if Pval.minute != Sval.minute:
            valor = 60 + Sval.second
            sr = str(valor) + sisec
        if len(sr)<4:
            sr = " " +  sr

    if flag ==1:
        sr = "    "
	# Construimos la oracion
    Phrases = Station + (5-len(Station))*" " + "  " + str(Pw) + " " + pstr + pr + "        " + sr  + "   " + str(Sw)

    return Phrases


# Paso a UTCDateTime de obspy desde el formato de hypoellipse
def hypoellipse_utcs(hypo_data):
    hypo_data = hypo_data.split("\n")[:-2]
    phases_values = []
    for phases in hypo_data:
        year = "20"+phases[9:11] + "-" + phases[11:13] + "-" + phases[13:15]+ "T"
        hour = phases[15:17] + ":" + phases[17:19]
        ptime = " "
        stime = " "

        if phases[7:8]   != ptime:
            ptime = obs.UTCDateTime(year + hour).timestamp + float(phases[20:22] + "." + phases[22:24])
        if phases[39:40] != stime:
            stime = obs.UTCDateTime(year + hour).timestamp + float(phases[32:34] + "." + phases[34:36])
        phases_values.append([ptime, stime])

    return np.array(phases_values)


#Para leer los resultados de hypoellipse del fichero .sum
def hypoellipse_reader(fichero):
    f = open(fichero)
    hypo = f.readlines()
    f.close()
    hypo = np.array(hypo, dtype = str)
    semiaxis4 = np.array([hypo[i][79:80] for i in range(len(hypo))])
    mask4 = semiaxis4 != "K"
    hypo=hypo[mask4]

    semiaxis1 = np.array([int(hypo[i][56:60]) for i in range(len(hypo))])
    semiaxis2 = np.array([int(hypo[i][60:63]) for i in range(len(hypo))])
    semiaxis3 = np.array([int(hypo[i][74:78]) for i in range(len(hypo))])

    #mask1 = semiaxis1!=9900
    #mask2 = semiaxis2!=9900
    #mask3 = semiaxis3!=9900
    #mask = mask1*mask2*mask3
    #hypo=hypo[mask]

    tiempos = np.array([hypo[i][:4] + "-" + hypo[i][4:6] + "-" + hypo[i][6:8] + "T" + hypo[i][8:10] + ":" + hypo[i][10:12] + ":" + hypo[i][12:14] + "." + hypo[i][14:16] for i in range(len(hypo))])
    longitudes = np.array([ -(float(hypo[i][24:26]) + float(hypo[i][27:29])/60 + float(hypo[i][29:31])/6000) if hypo[i][27:29] != '  ' else -(float(hypo[i][24:26]) +  float(hypo[i][29:31])/6000) for i in range(len(hypo))])
    latitudes = np.array([ float(hypo[i][16:18]) + float(hypo[i][19:21])/60 + float(hypo[i][21:23])/6000 if hypo[i][19:21] != '  ' else (float(hypo[i][16:18]) +  float(hypo[i][21:23])/6000) for i in range(len(hypo)) ])
    prof = np.array([ float(hypo[i][31:34]) + float(hypo[i][34:36])/100 if hypo[i][31:34] != '   ' else float(hypo[i][34:36])/100 for i in  range(len(hypo)) ])
    tiempos = np.array([obs.UTCDateTime(t.replace(" ", "0")).timestamp for t in tiempos])
    return tiempos, latitudes, longitudes, prof


# Para leer los resultados de hypoellipse en el Out
def hypoellipse_reader_out(fichero, return_fases = False ):
    """Funcion de lectura del .out de hypoellipse:
    In:
        fichero: Direccion del fichero ".out" de hypoellipse
    Out:
        parametros origen: array-like. tiempo origen, latitud, longitud, profundidad, rms, error horizontal, error vertical
    """

    with open(fichero) as f:
        hypo = np.array(f.readlines())


    # Cabeceras dentro del fichero
    eventos = np.arange(len(hypo))[hypo == '    date    origin      lat      long    depth    mag no d1 gap d  rms    avwt   se\n']+1
    errores = np.arange(len(hypo))[hypo == ' horizontal and vertical single variable standard deviations (68% - one degree of freedom; max 99 km)\n']
    
    # Parametros del origen del terremoto por hypoellipse
    tie, lat, lon, pro, rms, numero_fases, angle_gap = np.zeros((7, len(eventos))).astype(str)
    semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2 = np.zeros((5, len(eventos))).astype(str)

    for indice in range(eventos.shape[0]): 
        indice_fechas = eventos[indice]
        valores_origen  = hypo[indice_fechas]
        angle_gap[indice] = float(valores_origen[60:63])
        numero_fases[indice] = float(valores_origen[53:56])
        tie[indice] = valores_origen[1:5] + "/" + valores_origen[5:7] + "/" + valores_origen[7:9] + "T" + valores_origen[10:12] + ":" + valores_origen[12:14] + ":" + valores_origen[15:20]
        
        if valores_origen[33] == "w":
            indice_coordenada = -1

        if valores_origen[33] == "e":
            indice_coordenada = 1

        lat[indice] = float(hypo[indice_fechas+1][21:29])
        lon[indice] = indice_coordenada*float(hypo[indice_fechas+1][31:39])

        if valores_origen[40:46].replace(" ", "") != "":
            pro[indice] = float(valores_origen[40:46])
            
        rms[indice] = float(valores_origen[66:72])

        indice_errores = errores[indice] 
        semiaxis1[indice] = hypo[indice_errores+1][15:19]
        semiaxis2[indice] = hypo[indice_errores+1][40:44]
        semiaxis3[indice] = hypo[indice_errores+1][65:69]
        azimuth1[indice]  = hypo[indice_errores+2][14:20]
        azimuth2[indice]  = hypo[indice_errores+2][39:44]

    indice_fases = np.argwhere(hypo == '  stn c pha remk p p-sec s-sec resid  std-er   dist  azm ain    tc c vthk  ttob-ttcal-dlay-edly=resid rmk stn pha sources\n')[0]
    if return_fases:
        return tie, lat, lon, pro, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases, hypo, indice_fases
    else: 
        return tie, lat, lon, pro, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases


# Para leer los resultados de hypoellipse en el arc
def hypoellipse_to_catalog(fichero, nombre_salida = "/home/sysop/Documentos/echeyde/envio_echeyde/master_cluster_v1/prueba_fases_El_Hierro_VPVS_variable",
                            magniutdes = ""):
    with open(fichero) as f:
        valores_arc = f.readlines()

    localizacion_valores = valores_arc[0]
    valores_fases = valores_arc[1:]

    if len(localizacion_valores[31:36].replace(" ", "")) != 0:
        tie = localizacion_valores[0:4] + "/" + localizacion_valores[4:6] + "/" + localizacion_valores[6:8] + " " + localizacion_valores[8:10].replace(" ", "0") + ":"  + localizacion_valores[10:12].replace(" ", "0") + ":" + localizacion_valores[12:14].replace(" ", "0") + "." + localizacion_valores[14:16].replace(" ", "0")
        lat = str(round(float(localizacion_valores[16:18].replace(" ", "0")) + float(localizacion_valores[19:21].replace(" ", "0"))/60 + float(localizacion_valores[21:23].replace(" ", "0"))/6000,4))
        lon = str(round(float(localizacion_valores[24:26].replace(" ", "0")) + float(localizacion_valores[27:29].replace(" ", "0"))/60 + float(localizacion_valores[29:31].replace(" ", "0"))/6000,4))
        pro = str((float(localizacion_valores[31:36].replace(" ", "0"))/100))
        num_estaciones = localizacion_valores[38:41]
        rms = str(float(localizacion_valores[47:51])/100)


        sig_lat = "-"
        sig_lon = "-"

        if localizacion_valores[18:19] == "N":
            sig_lat = " "

        if localizacion_valores[26:27] == "E":
            sig_lon = " "


        lat = sig_lat + lat
        lon = sig_lon + lon
        mag = "    "

        if type(magniutdes) != str:
            mag = str(magniutdes)

        cabeza = "EVENTO\nFecha \t\t\t\t\t Latitud \t Longitud \t Pro \t No \t RMS \t Mag \n"
        cabeza += tie + "\t " + lat + "\t " + lon + "\t " + pro + "\t " + num_estaciones + "\t " + rms + "\t" + mag +"\n\n"

        phrase = ""
        cabecera = "sta	     fase   fecha                    error  calidad \n"
        phrase = cabecera

        for i in range(len(valores_fases))[:-2]:
            fase = valores_fases[i]
            Station = fase[:4]
            Pw = fase[7:8]
            Sw = fase[39:40]
            
            if fase[40:41] == "\n":
                continue

            fecha0 = obs.UTCDateTime("20" + fase[9:11] + "/" + fase[11:13] + "/" + fase[13:15] + " " + fase[15:17] + ":" + fase[17:19])

            if (Pw != "4") & (Pw != " "):
                fechaP = str(obs.UTCDateTime(fecha0.timestamp +float(fase[20:22]) +  float(fase[22:24])/100)).replace("T", " ")[:-5].replace("-", "/")
                errorP = float(fase[75:80])/100.

                if errorP>=0:
                    errorP = " " + str(errorP)

                else:
                    errorP = str(errorP)

                phrase += Station + "\t " + "P" + "\t " + fechaP + "\t\t" + errorP + "\t" + Pw + "\n"

            if (Sw != "4") & (Sw != " ") :
                fechaS = str(obs.UTCDateTime(fecha0.timestamp +float(fase[32:34]) +  float(fase[34:36])/100)).replace("T", " ")[:-5].replace("-", "/")
                errorS = float(fase[84:89])/100.

                if errorS>=0:
                    errorS = " " + str(errorS)

                else:
                    errorS = str(errorS)

                phrase += Station + "\t " + "S" + "\t " + fechaS + "\t\t" + errorS + "\t" + Sw + "\n"

        cabeza = cabeza + phrase + "\n"
        
        with open(nombre_salida, "a+") as f: 
            f.write( cabeza )


# Para arreglar los tiempos que salen de hypoellipse
def tiempos_fixer(tiempos):
    for i in range(len(tiempos)):
        tiempos[i]=tiempos[i].replace(" ", "0")
        try:
            tiempos[i] = tiempos[i].replace("Z","")
        except:
            None
        if float(tiempos[i].split(":")[-1].replace(" ","0"))>59.99:
            tiem = tiempos[i].split(":")
            tiem[-2] = str(int(tiem[-2])+1)
            tiempos[i] = tiem[0]+":"+tiem[1]+ ":"+"00.00"
    return tiempos


def hypoellipse_vpvs_change(name, fichero_modelo, fichero_vp_vs):
    fecha, hora, ratio, poisson = np.genfromtxt(fichero_vp_vs, dtype = str).T
    time_stamps = np.array([obs.UTCDateTime(fecha[i] + " " + hora[i]).timestamp for i in range(len(hora))])
    tiempo = float(name.replace("_tmp", ""))
    f2 = interp1d(time_stamps, ratio, kind = "cubic")   

    if min(time_stamps)>tiempo:
        vel = ratio[0]

    if (min(time_stamps) < tiempo) & (max(time_stamps)> tiempo):
        vel = f2(tiempo).round(2)

    if max(time_stamps) < tiempo : 
        vel = ratio[-1]
        
    with open(fichero_modelo, "r") as fic, open(name + ".prm", "w") as fic2:
        lines = fic.readlines()
        phrase = ""
        for line in lines:
            if line[:5] == "VELOC":
                line = line[:-5] + str(vel)

            phrase += line + "\n"
        fic2.writelines(phrase)
    

# Localizador con hypoellipse
def hypoellipse_locator(hypo_data, name_file = "localizacion", hypoellipse_route = "/home/sysop/Documentos/echeyde/envio_echeyde/programas_localizacion/hypoellipse3", \
                            fichero_vp_vs = None, hypoin_file = 'hypo_hierro.in', hypoctl_file = 'hypoc_hierro.ctl', remove = True,  nombre_salida = "prueba", \
                            return_values = True, write_catalog_file = True,  magnitudes = ""):

    """
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
    name = str(name_file) + "_tmp"
    b = os.popen("pwd").readline()[:-1]

    #Vamos al directorio de hypoellipse
    os.chdir(hypoellipse_route)
    with open(name + ".pha", "w") as f: 
        f.write(hypo_data)

    #Reescribimos el hypo.in para nuestro fichero
    with open(hypoin_file, "r") as fd1, open(name + hypoin_file, "w") as fd1_0: 
        linea1 = fd1.readlines()
        linea1[0] = name + hypoctl_file + '\n'
        linea1[1] = "yes\n"
        linea1[2] = name + '\n'
        fd1_0.writelines(linea1)

    # Escribimos el hypoctl para nuestro fichero
    with open(hypoctl_file, "r") as fd2, open(name + hypoctl_file, "w") as fd2_0:
        linea2 = fd2.readlines()
        linea2[-2] = 'jump ' + name + ".pha\n"
        fichero_modelo = linea2[9].replace('jump ', "").replace("\n", "")

        if fichero_vp_vs:
            linea2[9] = 'jump ' + name + ".prm\n"
        
        fd2_0.writelines(linea2)

    if fichero_vp_vs:
        hypoellipse_vpvs_change(name, fichero_modelo, fichero_vp_vs)


    # Localizamos el evento
    os.system("./hypoel < " + name + hypoin_file + " >> tmp.txt")

    
    if write_catalog_file:
        
        hypoellipse_to_catalog(name + ".arc", nombre_salida, magnitudes)

    if remove:
        os.system("rm *"+ name + "*")

    os.chdir(b)
    if return_values: 
        #Extraemos los valores que queremos
        tiempos, latitudes, longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases = hypoellipse_reader_out(name+".out")
        return tiempos, latitudes, longitudes, prof, rms, semiaxis1, semiaxis2, semiaxis3, azimuth1, azimuth2, angle_gap, numero_fases

