#from ensemblrest import EnsemblRest
import urllib2
import simplejson as json
import sys
from Bio import SeqIO
import cStringIO as StringIO

def cds_quality_check(seq):
	##CHECK if sequence is inframe, has start & stop codons
	##USED TO CHECK STOP CODONS ##NOT IMPLEMENTED YET
	stop_codons = ["TGA", "TAG", "TAA"]
	#is_stop_codon = lambda x : any(x[i:i+3] in stop_codons for i in range(0,len(x),3))
	inframe_check=start_check=stop_check=False ##Assume that all checks will be false 
	stop_codons
	seq_as_file=StringIO.StringIO(seq)
	fasta = SeqIO.parse(seq_as_file, "fasta")
	for i , seq_record 	in enumerate(fasta):
		seq_len=len(seq_record.seq)
		inframe_check=seq_len%3==0
		start_check=seq_record.seq[0:3]==str("ATG")
		stop_check=any(x in seq_record.seq[seq_len-3:seq_len] for x in stop_codons)
	if(all([inframe_check,start_check,stop_check])):
		print seq
		sys.exit(0)
	else:
		sys.exit(2)

##ENTRYPOINT

gene_id=sys.argv[1]

#ensRest = EnsemblRest()

try:
	##http://rest.ensembl.org/lookup/id/Bra032705?expand=1;content-type=application/json
	base_lookup="http://rest.ensembl.org/lookup/id/"+gene_id+"?expand=1;content-type=application/json"
	#gene_info=ensRest.getLookupById(id=gene_id)
	#print gene_info
	output=urllib2.urlopen(base_lookup).read()
	lookup_out=json.loads(output)
	transcript_id=(json.dumps(lookup_out["Transcript"][0]['id'])).strip('"') ##GET ONLY MAIN PROTEIN CODING SEQUENCE
	base_seq='http://rest.ensembl.org/sequence/id/'+transcript_id+'?type=cds;object_type=transcript;content-type=text/x-fasta'
	cds=urllib2.urlopen(base_seq).read()
	#gene_seq=ensRest.getSequenceById(id=gene_info['id'],feature=["protein_coding","cds","exon"])
	#print gene_seq
	pass
except:
	sys.exit(2) ##CANT FIND GENE

#print ">",gene_info['species'],"|",gene_info['id']
#print gene_seq['seq']
cds_quality_check(cds)

sys.exit(2) ##IF didnt exit after quality check then there's some error in sequence, so exit 2