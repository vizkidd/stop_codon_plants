import csv
import sys
import re

def findWholeWord(w):
    return re.compile(r'\b({0})\b'.format(w), flags=re.IGNORECASE).search

##ENTRYPOINT

l2s_file =  sys.argv[1]
og_file=sys.argv[2] #UNUSED
og2_file=sys.argv[3] #UNUSED
tax_id=sys.argv[4] 

# using the tax_id, extract org_id from levels2species. (extracting plant organisms)
#*if third column contains the tax_id then we store the 2nd column(org_id) of levels2species
f=open(l2s_file, 'r')
csv_reader = csv.reader(f, delimiter='\t')
for row in csv_reader:
	if(findWholeWord(tax_id)(row[3]) is not None):#(tax_id in row[3]):
		print(str(row[1]))
f.close()

sys.exit(0)
