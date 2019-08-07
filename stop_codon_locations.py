import sys
from Bio import SeqIO
import pickle

def stop_codon_locations(seq):
	for stop_codon in stop_codons:
		stop_pos=str(str(seq.seq).find(stop_codon))
		ele=seq.id+"\t"+stop_codon+"\t"+stop_pos
		#print(ele)
		if(not str(stop_pos) in stop_dict):
			stop_dict[stop_pos]=[]
		stop_dict[stop_pos].append(seq.id)

def unique(list1): 
    # intilize a null list 
    unique_list = []  
    # traverse for all elements 
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    # print list 
    for x in unique_list: 
        print x, 

##ENTRYPOINT
##GIVE fasta FILE & mol_type AS PARAMETERS

seq_file=sys.argv[1]
mol_type=sys.argv[2]
out_file=sys.argv[3]

stop_dict =	{}

if(mol_type=="DNA"):
	stop_codons = ["TGA", "TAG", "TAA"]
elif(mol_type=="RNA"):
	stop_codons = ["UGA", "UAG", "UAA"]
elif(mol_type=="AA"):
	stop_codons = ["*"]
else:
	sys.exit(-1)

#print("Stop Codons:")
#print(stop_codons)

rec_count=0
fasta = SeqIO.parse(seq_file, "fasta")
for record in fasta:
	stop_codon_locations(record)
	rec_count=rec_count+1
      
pickle_out = open(out_file,"wb")
pickle.dump(stop_dict, pickle_out)
pickle_out.close()

		