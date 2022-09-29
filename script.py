from Bio import SeqIO
from Bio import Seq

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

#Compara los largos de las secuencias dentro de una lista y se queda con la más larga (ESTA!)	
def winner(lista):
	l=0
	for i in range(len(lista)):
		if len(lista[i])>l:
			messi=lista[i]
			l=len(lista[i])
	return messi

#Importo archivo .gb (es ADN) 
DNA=SeqIO.read("NTRK1.gb","genbank")
#COn reverse_complement() genero la cadena de ARN complementaria
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


