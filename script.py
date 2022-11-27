from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Align.Applications import MuscleCommandline
import os

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

def callblast(archivo):
    record = SeqIO.read(archivo, format="fasta")
    result_handle = NCBIWWW.qblast("blastp","nr", record.format('fasta'),descriptions=100,alignments=100) 
    blast=open('blast.xml','w')
    blast.write(result_handle.read())
    blast.close()
    return 

def MSA(archivo):
	#Creo la variable con las alineaciones
	secuencias_aln=MuscleCommandline(input=archivo)

	#Genero un archivo FASTA del output
	seq_aln=open('seq_aln.fasta','w')
	for secuencias_aln in SeqIO.parse(archivo,"fasta"):
		seq_aln.write(">")
		seq_aln.write(str(secuencias_aln.id))
		seq_aln.write(" ")
		seq_aln.write(str(secuencias_aln.description))
		seq_aln.write("\n")
		seq_aln.write(str(secuencias_aln.seq))
		seq_aln.write("\n")
	seq_aln.close()
	return 

#Importo archivo .gb (es ADN) 
DNA=SeqIO.read("NTRK1.gb","genbank")
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
SeqIO.write(DNA,"protein.fasta","fasta")

#Blast
#callblast("protein.fasta")

#MSA 

MSA("secuencias.fasta")

#Uso EMBOSS (función getorf) para calcular los ORFs
os.system("getorf NTRK1.gb ORFs.fasta -find=1")

#Comparo con prosite. Para eso uso las funciones de EMBOSS prosextract y patmatmotifs
os.system("prosextract ./")
os.system("patmatmotifs protein.fasta dominios")

#BUSCADOR 

#Ingreso la búsqueda
a_buscar=input("Ingrese búsqueda:")

#Abro el archivo XML donde voy a realizar la búsqueda
archivo=open("blast.xml","r")
archivo_guardado= NCBIXML.parse(archivo)
item=next(archivo_guardado)

#En caso de que se haya ya ejecutado el script borro el archivo Ex4.txt para que se cree uno nuevo
if path.exists('Ex4.txt'):
    remove('Ex4.txt')

#Busco (no diferencia entre mayúsculas y minúsculas) y guardo en el archivo Ex4
for alignment in item.alignments:
          for hsp in alignment.hsps:
              if alignment.hit_def.find(a_buscar.casefold()) != -1:
                  Ex4=open('Ex4.txt','a')
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
                  Ex4.close()
                  """
                  Entrez.email = 'agortiz@itba.edu.ar'
                  handle = Entrez.efetch(db="protein", id=alignment.hit_id, rettype="fasta")
                  print(handle.read())
                  """

