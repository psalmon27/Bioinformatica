#MÓDULOS y LIBRERÍAS


from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import os
from os import remove
from os import path
import configparser

#FUNCIONES


#Recibe una secuencia de ARN y genera los dos ORF que faltan (2 y 3) y los traduce a proteina. RECORDAR que el archivo que importamos nos da la info para el primer ORF.
def ORFfinder(sec):

	orf2=sec[1:]
	orf3=sec[2:]
	
	orf2=orf2.translate()
	orf3=orf3.translate()
	
	return orf2,orf3

#Recibe los ORFs traducidos a aa. Busca donde hay una M (inicio) y donde hay un * (stop). Acto seguido, los splitea y guarda los splits en una lista.
def splitter(sec):

	ORFs=[]
	start=sec.find("M")
	end=0
	while start>=0 and end>=0:
		sec=sec[start:]
		end=sec.find("*")
		split=sec[0:end]
		ORFs.append(split)
		sec=sec[end+1:]
		start=sec.find("M")
	return ORFs

#Compara los largos de las secuencias dentro de una lista y se queda con la más larga	
def winner(lista):
	l=0
	for i in range(len(lista)):
		if len(lista[i])>l:
			messi=lista[i]
			l=len(lista[i])
	return messi

#Ejecuta blast online

def callblast(query,tipo,bd,resultados,matriz,e_thres):
    record = SeqIO.read(query, format="fasta")
    result_handle = NCBIWWW.qblast(tipo,bd, record.format('fasta'),descriptions=resultados,alignments=resultados,hitlist_size=resultados,matrix_name=matriz,expect=e_thres) 
    blast=open('blast.xml','w')
    blast.write(result_handle.read())
    blast.close()
    return 


#Buscador de patterns

def buscador(pattern,query):

	#Abro el archivo XML donde voy a realizar la búsqueda
	archivo=open(query,"r")
	archivo_guardado= NCBIXML.parse(archivo)
	item=next(archivo_guardado)

	#En caso de que se haya ya ejecutado el script borro el archivo Ex4.txt para que se cree uno nuevo
	if path.exists('hits.txt'):
		   remove('hits.txt')

		#Busco (no diferencia entre mayúsculas y minúsculas) y guardo en el archivo Ex4
	for alignment in item.alignments:
		for hsp in alignment.hsps:
			if alignment.hit_def.casefold().find(pattern.casefold()) != -1:
				Ex4=open('hits.txt','a')
				Ex4.write('****Alignment****')
				Ex4.write("\n")
				Ex4.write('sequence: ')
				Ex4.write(alignment.title)
				Ex4.write("\n")
				Ex4.write('length: ')
				Ex4.write(str(alignment.length))
				Ex4.write("\n")
				Ex4.write('score: ')
				Ex4.write(str(hsp.score))
				Ex4.write("\n\n")
				
				Ex4.write('FASTA: ')
				Entrez.email = 'agortiz@itba.edu.ar'
				handle = Entrez.efetch(db="protein", id=alignment.hit_id, rettype="fasta")
				Ex4.write(str(handle.read()))
				Ex4.write("\n\n")
				Ex4.close()
				
	return


#CÓDIGO PRINCIPAL


#Importo archivo de configuración
configuracion=configparser.ConfigParser()
configuracion.read('entradas.cfg')

#1: Búsqueda y selección de ORFs

#Importo archivo .gb (es ADN) 
NTRK1=configuracion['Búsqueda de ORFs']['archivo']+".gb"
DNA=SeqIO.read(NTRK1,"genbank")
#Con reverse_complement() genero la cadena de ARN complementaria
rna_comp=DNA.reverse_complement().seq
#Transcribo el ADN (cadena lider) a ARN
DNA.seq=DNA.seq.transcribe()

#Los leader son de la secuencia original y los comp son de la cadena complementaria (en total son 6)
#Para leader2/3 y comp2/3 es necesario usar la función ORFfinder
leader1=DNA.seq.translate()
comp1=rna_comp.translate()
leader2,leader3=ORFfinder(DNA.seq)
comp2,comp3=ORFfinder(rna_comp)

#Genero listas con los splits de todas las cadenas leader y comp
frames1=splitter(leader1)
frames2=splitter(comp1)
frames3=splitter(leader2)
frames4=splitter(comp2)
frames5=splitter(leader3)
frames6=splitter(comp3)

#Buscamos la secuencia con mayor longitud dentro de cada frame
op1=winner(frames1)
op2=winner(frames2)
op3=winner(frames3)
op4=winner(frames4)
op5=winner(frames5)
op6=winner(frames6)

#Buscamos la secuencia con mayor longitud de todas las op
diego=winner([op1,op2,op3,op4,op5,op6])

#Modificamos el atributo .seq de DNA por la proteína diego
DNA.seq=diego

#Guardo ambas secuencias en formato fasta
SeqIO.write(DNA,"bestORF.fasta","fasta")

#Uso EMBOSS (función getorf) para calcular los ORFs (ejercicio 5)
os.system("getorf "+NTRK1+" ORFs.fasta -find=1")

#2. Blast (la función está comentada porque tarda aprox 2 minutos en ejecutarse)

query=configuracion['Blast']['archivo']+".fasta"
tipo=configuracion['Blast']['tipo']
bd=configuracion['Blast']['bd']
resultados=int(configuracion['Blast']['resultados'])
matriz=(configuracion['Blast']['matriz'])
e_thres=float(configuracion['Blast']['E umbral'])
callblast(query,tipo,bd,resultados,matriz,e_thres)

#3. MSA 

input_msa=configuracion['MSA']['archivo']+".fasta"
os.system("muscle -in "+input_msa+" -out seq_aln.fasta")

#4. Dominios

proteina=configuracion['Dominios']['archivo']+".fasta"
#Comparo la proteína con prosite. Para eso uso las funciones de EMBOSS prosextract y patmatmotifs
os.system("prosextract ./")
os.system("patmatmotifs "+proteina+" dominios")

#5. Búsqueda de pattern

#Ejecuto buscador de patterns
pattern = configuracion['Búsqueda de pattern']['pattern']
query = configuracion['Búsqueda de pattern']['archivo']+".xml"
buscador(pattern,query)
