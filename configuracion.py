import configparser
configuracion = configparser.ConfigParser()

configuracion['Búsqueda de ORFs'] = {'archivo' : 'NTRK1'}
configuracion['Blast'] = {'archivo' : 'protein', 'Tipo' : 'blastp', 'BD' : 'nr', 'Resultados' : '100'}
configuracion['MSA'] = {'archivo' : 'secuencias'}
configuracion['Búsqueda de pattern'] = {'archivo' : 'blast', 'pattern' : 'sapiens'}
configuracion ['Dominios'] = {'archivo' : 'protein'}


with open('entradas.cfg', 'w') as archivoconfig:
    configuracion.write(archivoconfig)

