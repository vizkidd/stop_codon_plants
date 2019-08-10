from intermine import registry
from intermine.webservice import Service
import intermine.results
import intermine.errors
import sys
import time

#registry.getInfo("phytomine")
#registry.getData("phytomine")

def cds_quality_check(seq,id_str):
	##CHECK if sequence is inframe, has start & stop codons
	##USED TO CHECK STOP CODONS ##NOT IMPLEMENTED YET
	stop_codons = ["TGA", "TAG", "TAA"]
	#is_stop_codon = lambda x : any(x[i:i+3] in stop_codons for i in range(0,len(x),3))
	seq_len=len(seq)
	inframe_check=seq_len%3==0
	start_check=seq[0:3]==str("ATG")
	stop_check=any(x in seq[seq_len-3:seq_len] for x in stop_codons)
	if(all([inframe_check,start_check,stop_check])):
		print id_str
		print seq
		sys.exit(0)
	else:
		sys.exit(2)

##https://github.com/k821209/pipelines/wiki/python.phytomine
def get_cds(name,id_str): # name should be transcript name
    template = service.get_template('Transcript_CDS_sequence')
    rows = template.rows(
        A = {"op": "=", "value": name}
        )
    for row in rows:
        cds_quality_check(row["CDSs.sequence.residues"],id_str)

def get_protein(id,id_str):
	template = service.get_template('Protein_protein_sequence')
	rows = template.rows(
		A = {"op": "=", "value": id}
		)
	for row in rows:
		print id_str
		print row["sequence.residues"]

#get_protein(row["proteins.primaryIdentifier"],id_str)
#"proteins.primaryIdentifier", "proteins.length

def query_db(genes_found):
		for gene in genes_found:
			#print gene
			query = service.new_query("Transcript")
			query.add_view(
				"CDSs.transcript.CDSs.primaryIdentifier", "name", "gene.name",
    			"protein.name", "organism.name"
			)
			query.add_constraint("gene.name", "=", gene, code = "A")
			query.add_constraint("primaryTranscript", "=", "true", code = "C")
			query.add_constraint("gene.secondaryIdentifier", "CONTAINS", gene, code = "B")
			query.set_logic("(A or B) and C")
			for row in query.rows():
				id_str=">"+row["CDSs.transcript.CDSs.primaryIdentifier"]+" | "+row["gene.name"]+" "+row["organism.name"]+" "+ row["name"]+" " +row["protein.name"]
				get_cds(row["name"],id_str)
				print "\n"

##ENTRYPOINT
##GIVE gene accession as PARAM

gene_acc=sys.argv[1]
#mol_type=sys.argv[2]

service = Service("https://phytozome.jgi.doe.gov/phytomine/service")
query = service.new_query("Transcript")
query.add_view("gene.name")
query.add_constraint("gene", "LOOKUP", gene_acc, code = "B")
query.add_constraint("primaryTranscript", "=", "true", code = "A")

genes_found=[]
try:
	for row in query.rows():
		#print row
		genes_found.append(row["gene.name"])
except intermine.errors.WebserviceError as e:
	time.sleep(3)
	for row in query.rows():
		#print row
		genes_found.append(row["gene.name"])

#print genes_found


if not genes_found:
	genes_found.append(gene_acc) ##gene_acc is probably a phytozome ID
	#sys.exit(2) ##EXIT BECASE WE CANT FIND GENE IN PHYTOZOME

try:
	query_db(genes_found)
except intermine.errors.WebserviceError as e:
	time.sleep(3)
	query_db(genes_found)

sys.exit(2) ##IF didnt exit after quality check then there's some error in sequence, so exit 2

#query.add_constraint("gene.secondaryIdentifier", "CONTAINS", "37374315", code = "B")