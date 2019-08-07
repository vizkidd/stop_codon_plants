import sys
from Bio import SeqIO


def trim(seq,len):
	print(">"+str(seq.id)[:trim_length].upper())
	print(seq.seq)
##ENTRYPOINT
##GIVE fasta FILE & trim length AS PARAMETERS

seq_file=sys.argv[1]
trim_length=int(sys.argv[2])

fasta = SeqIO.parse(seq_file, "fasta")
for record in fasta:
	trim(record,trim_length)
