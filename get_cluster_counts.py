import csv
import sys
from itertools import chain

file=sys.argv[1]
lower_limit=int(sys.argv[2])
#cluster_file=sys.argv[3]
fieldnames=['cluster', 'genes']
#cluster_list=csv.reader(open(cluster_file))
#clusters=list(chain.from_iterable(cluster_list))
#print clusters
with open(file) as csvfile:
	reader = csv.DictReader(csvfile,fieldnames=fieldnames,delimiter='\t')
	#IT is iterating rowwise so we need to store clusterwise data
	temp_list=[]
	curr_clust=reader.next()['cluster']
	#print(type(curr_clust))
	for row in reader:
		prev_clust=curr_clust
		curr_clust=row['cluster']
		if(prev_clust==curr_clust):
			#Add data to a temp dict
			temp_list.append(row['genes'])
		else:
			#if(str(prev_clust) in clusters):
				##prcoess data since we've moved away from the cluster
			cnt=len(temp_list)
			if(cnt>=lower_limit):
				print prev_clust+"\t"+str(cnt)
			temp_list=[]
			temp_list.append(row['genes'])
