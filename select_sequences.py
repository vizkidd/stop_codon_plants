import pickle
import sys
from Bio import SeqIO
import statistics
from operator import itemgetter
import random

seq_file=sys.argv[1]
name_file=sys.argv[2]
out_file=sys.argv[3]

pickle_in = open(name_file,"rb")
name_dict = pickle.load(pickle_in)

seq_dict={}
recs = SeqIO.index(seq_file, "fasta")
fasta = SeqIO.parse(seq_file, "fasta")
len_list=[]
for record in fasta:
	#if(record.name.split(':')[0] in name_dict):
	len_list.append(len(record.seq))
	org_id=record.name.split(':')[0]
	if(org_id not in seq_dict):
		seq_dict[org_id]=[]
	seq_dict[org_id].append(record)
		
avg_len=statistics.mean(len_list)
#print len_list
print "Mean Len:"+str(avg_len)

##Calculate deviation
dev_dict={}
for key,val in seq_dict.items():
	#if(len(val)>1):
	dev_dict[key]=[]
	for ele in val:
		dev_dict[key].append((avg_len-len(ele.seq))**2)

chosen_seqs={}
##FIND the index of lengths with min deviation
for key,val in dev_dict.items():
	if(key in name_dict):
		if(len(val)>1):
			##if all elements are same then choose randomly
			if(all(elem == val[0] for elem in val)):
				chosen_seqs[key]=seq_dict[key][random.randint(0,len(val))-1]
			##else choose sequence with minimum deviation
			else:
				chosen_seqs[key]=seq_dict[key][min(enumerate(val), key=itemgetter(1))[0]]
		else:
			chosen_seqs[key]=seq_dict[key][0]
		chosen_seqs[key].id=name_dict[key]
		chosen_seqs[key].description=""

#print list(chosen_seqs)
SeqIO.write([ v for v in chosen_seqs.values() ],out_file,"fasta")
