import configparser
configuracion = configparser.ConfigParser()

configuracion['Búsqueda de ORFs (archivo en formato Genbank)'] = {'Nombre del archivo' : 'NTRK1'}
configuracion['Blast (archivo en formato FASTA)'] = {'Nombre del archivo' : 'protein', 'Tipo' : 'blastp', 'Base de datos' : 'nr', 'Cantidad de resultados' : '100'}
configuracion['MSA (archivo en formato FASTM)'] = {'Nombre del archivo' : 'secuencias'}
configuracion['Búsqueda de pattern (archivo en formato XML)'] = {'Nombre del archivo' : 'blast', 'pattern' : 'sapiens'}
configuracion ['Dominios (archivo en formato FASTA)'] = {'Archivo' : 'protein'}


with open('entradas.cfg', 'w') as archivoconfig:
    configuracion.write(archivoconfig)

