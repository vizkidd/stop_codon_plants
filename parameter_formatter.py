import sys

pars=[]
pars=sys.argv[1:len(sys.argv)]

print(" ".join(str(x) for x in pars))
