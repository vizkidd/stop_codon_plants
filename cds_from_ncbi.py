import urllib2
import sys, getopt
#import ElementTree
import xml.etree.ElementTree as ET
from Bio import SeqIO
import cStringIO as StringIO
from Bio import Entrez
import time
import re
from lxml import objectify

MAX_RETRIES=4
RETRY_COUNT=0

def findWholeWord(w):
    return re.compile(r'{0}'.format(w), flags=re.IGNORECASE).search

def get_cds(gid,species):
	global MAX_RETRIES, RETRY_COUNT
	#web_env1 = obj.find('WebEnv').text
	#key_1 = obj.find('QueryKey').text
	if(RETRY_COUNT>=MAX_RETRIES):
		return
	base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	url = base + "elink.fcgi?dbfrom="+"gene"+"&db="+db2+"&id="+gid
	#url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom="+db1+"&dbto="+db2+"?LinkName=gene_nuccore&from_uid="+gid
	output = urllib2.urlopen(url).read()
	#print output

	obj = objectify.fromstring(output)
	try:
		lnk_len=len(obj.LinkSet.LinkSetDb.Link)
		for lnk in obj.LinkSet.LinkSetDb.Link:
			##get id elements which are two nodes deep
			url = base + "efetch.fcgi?db="+db2+"&id="+str(lnk.Id);
			url = url + "&rettype=fasta_cds_na&retmode=text";
			cds = urllib2.urlopen(url).read()
			seq_as_file=StringIO.StringIO(cds)
			fasta = SeqIO.parse(seq_as_file, "fasta")
			if(len(list(fasta))==1):
				cds_quality_check(cds,species)
			time.sleep(0.3)
	except AttributeError as e:
		sys.exit(2)
	#RETRY_COUNT+=1
	#time.sleep(1)
	#if(mol=="protein"):
	#	get_cds(gid,species,"gene")
	#else:
	#	get_cds(gid,species,"protein")

def cds_quality_check(seq,species):
	##CHECK if sequence is inframe, has start & stop codons
	##USED TO CHECK STOP CODONS ##NOT IMPLEMENTED YET
	stop_codons = ["TGA", "TAG", "TAA"]
	#is_stop_codon = lambda x : any(x[i:i+3] in stop_codons for i in range(0,len(x),3))
	inframe_check=start_check=stop_check=False ##Assume that all checks will be false 
	#stop_codons
	seq_as_file=StringIO.StringIO(seq)
	fasta = SeqIO.parse(seq_as_file, "fasta")
	for i , seq_record in enumerate(fasta):
		name_list=str(seq_record.name).split("_")
		if(name_check):
			if(name not in seq_record.name):#if(findWholeWord(name)(seq_record.name) is None):
				continue
		seq_len=len(seq_record.seq)
		#print(seq_len)
		inframe_check=seq_len%3==0
		#start_check=seq_record.seq[0:3]==str("ATG")
		stop_check=any(x in seq_record.seq[seq_len-3:seq_len] for x in stop_codons)
		#print str(inframe_check)+" "+str(stop_check)
		if(all([inframe_check,stop_check]) and len(str(seq_record.seq))>1): #start_check
			print ">"+species+"_"+name_list[len(name_list)-2]+"_"+name_list[len(name_list)-1]+" "+seq_record.name
			print seq_record.seq
			sys.exit(0)

##ENTRYPOINT

##PASS gene accession AS PARAM
db2 = 'nuccore'

if(len(sys.argv)!=3):
	sys.exit(2)
id_list = sys.argv[1] #'At1g23480'
species=sys.argv[2]
#db1=sys.argv[3]
#print(id_list, species, db1)
name_check=(not id_list.isdigit())
name=id_list
linkname ="gene_nucleotide"
#print name_check
#print name

base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
url = base + "esearch.fcgi?db="+"gene"+"&term="+id_list+"&usehistory=y&retmax=1"

output=urllib2.urlopen(url).read()
#print output

#print "------------------------------"

obj = ET.fromstring(output)

first_id=None

id_list=[]

if obj.find('IdList') is None: ##CANT FIND GENE
	sys.exit(2)

for child in obj.find('IdList'):
	id_list.append(child.text)

if len(id_list) == 0:  ##CANT FIND GENE
	sys.exit(2)
#print(id_list)
try:
        get_cds(id_list[0],species)
except urllib2.HTTPError as e:
	RETRY_COUNT+=1
        time.sleep(3)
        get_cds(id_list[0],species)
##IF didnt exit after quality check then there's some error in sequence, so exit 2
sys.exit(2)
##OLD CODE
#for child in obj.iter('Id'):
#    first_id=child.text
#    break

#if first_id is None: 
#	sys.exit(2)

#print is_stop_codon(sequence)
#print is_stop_codon(sequence2)

##DOESNT WORK BECAUSE BIOPYHTON DOESNT GET CDS
##Biopython E-fetch
#Entrez.email = "vishveshkarthik@gmail.com" #required
#seq_record = SeqIO.read(Entrez.efetch(db=db2, rettype="gb", retmode="text", id=ext_id), "gb")

#print seq_record

#gb_record = SeqIO.read(genbank, "genbank")

#print gb_record.name
#print len(gb_record.features)

