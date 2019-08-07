import sys
from Bio import SeqIO
from Bio import Seq
import pickle
import collections

def cut(seq):
	seq_len=len(seq.seq)
	up_seq=str(seq.seq).upper()
	stp=up_seq[seq_len-3:seq_len]
	stop_check=(stp in stop_codons)
	#print(stop_check)
	cdns=[stp == i for i in stop_codons]
	if(stop_check):
		for cdn in range(3): #which(stop_check)
			if(cdns[cdn]):
				return cdn+1
	else:
		return "NA"	

##ENTRYPOINT
##GIVE fasta FILE & trim length AS PARAMETERS

seq_file=sys.argv[1]
cluster=sys.argv[2]
new_file=sys.argv[3]
name_file=sys.argv[4]
min_nas=sys.argv[5]

stop_codons = ["TAG", "TGA", "TAA"]
row={}
recs=[]

pickle_in = open(name_file,"rb")
name_dict = pickle.load(pickle_in)

#row["name"]=str(cluster)
for key, value in name_dict.items():
	row[value]="NA"

fasta = SeqIO.parse(seq_file, "fasta")
for record in fasta:
	if(record.name in row):
		row[record.name]=str(cut(record))
	else:
		continue
	record.seq=record.seq[:len(record.seq)-3].upper()
	for cdn in range(0,len(record.seq)-3,3):
		trimer=str(record.seq[cdn:cdn+3])
		for i in stop_codons:
			if(i==trimer):
				record.seq=Seq.Seq("".join((str(record.seq[:cdn]),"NNN",str(record.seq[cdn+3:]))))
	recs.append(record)
	SeqIO.write(recs,new_file,"fasta")

row_ord=collections.OrderedDict(sorted(row.items()))
#        row_ord[key]=value

#print row_ord
#print(str(cluster)+" "+" ".join(str(x) for x in row_ord.keys()))
if(row_ord.values().count("NA")<=int(min_nas)):
	print(str(cluster)+" "+" ".join(str(x) for x in row_ord.values()))
#print(*row, sep=' ')
