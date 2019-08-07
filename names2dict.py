import csv
import pickle
import sys
from itertools import chain
import collections

names_file=sys.argv[1]
fieldnames = ['org', 'name']
names=csv.reader(open(names_file),delimiter='\t')#, fieldnames=fieldnames)
names_dict={}
for row in names:
	names_dict[row[0]]=row[1]
#print names_dict

names_list=collections.OrderedDict()
for key,value in sorted(names_dict.items(),key=lambda item:item[1]):
	names_list[key]=value.upper()

#print names_list
pickle_out = open(names_file+".dict","wb")
pickle.dump(names_list, pickle_out)
pickle_out.close()

print(" ".join(str(x) for x in names_list.values()))
