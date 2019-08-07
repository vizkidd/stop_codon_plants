import sys
from Bio import SeqIO
import pickle


seq_file=sys.argv[1]
in_file=sys.argv[2]
out_dir=sys.argv[3]
lower_limit=int(sys.argv[4])

pickle_in = open(in_file,"rb")
stop_dict = pickle.load(pickle_in)
file_list=[]

#print(out_dir)
recs = SeqIO.index(seq_file, "fasta")
#print(recs)
##Iterate through dict and group the sequences into clusters based on stop codons
for key,ls in stop_dict.items():
	clust_list=[]
	file_list.append(key)
	for seq_id in ls:
		if(str(seq_id) in recs):
			clust_list.append(recs[seq_id])
	#print(key)
	#print(clust_list)
	if(len(clust_list)>lower_limit):
		print(key)
		SeqIO.write(clust_list,out_dir+"/"+key+".fasta","fasta")

#pickle_out = open(out_dir+"/filelist.list","wb")
#pickle.dump(file_list, pickle_out)
#pickle_out.close()
