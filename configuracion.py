import configparser
configuracion = configparser.ConfigParser()

configuracion['Búsqueda de ORFs'] = {'archivo' : 'NTRK1'}
configuracion['Blast'] = {'archivo' : 'bestORF', 'Tipo' : 'blastp', 'BD' : 'nr', 'Resultados' : '100', 'Matriz' : 'BLOSUM62', 'E umbral' : '0.05'}
configuracion['MSA'] = {'archivo' : 'secuencias'}
configuracion['Búsqueda de pattern'] = {'archivo' : 'blast', 'pattern' : 'sapiens'}
configuracion ['Dominios'] = {'archivo' : 'bestORF'}


with open('entradas.cfg', 'w') as archivoconfig:
    configuracion.write(archivoconfig)
