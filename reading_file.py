import pandas as pn
a=pn.read_csv('GC_quant.bed',delimiter="\t",header=0)

f=a[(a.ix[:,7]>=0.35) & (a.ix[:,7]<=0.75) & (a.ix[:,14]==20)]

f.to_csv('remaining.bed',sep='\t',index=False)
