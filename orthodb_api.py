import demjson
import urllib2
import xml.etree.ElementTree as ET
from Bio import SeqIO
import cStringIO as StringIO
from Bio import Entrez
import sys
import re
import csv

def fetch_data(gid):
	url = "http://www.orthodb.org/ogdetails?id="+gid;
	output = urllib2.urlopen(url).read()
	js=demjson.decode(output)
	#print js
	#ncbi_id=js["data"]["entrez"][0]["id"]
	#print ncbi_id
	try:
		id_val=""
		if("entrezprotein" in js["data"]):
			id_val=js["data"]["entrezprotein"][0]["id"] ##protein
		if("xrefs" in js["data"]):
			for keys in js["data"]["xrefs"]:
				if(keys["type"]=="GeneID"):
					id_val=keys["id"] ##gene
		if("entrez" in js["data"]):
			id_val=js["data"]["entrez"][0]["id"] #protein
		if(id_val!=""):
			print(id_val)
		else:
			sys.exit(2)
		sys.exit(0)
	except KeyError as e:
		print("-1"+"\t"+"-1")
		sys.exit(2)


##ENTRYPOINT

gid=sys.argv[1]

fetch_data(gid)
sys.exit(2)
