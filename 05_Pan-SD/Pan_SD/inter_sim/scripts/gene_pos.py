import json 
import re
import pandas as pd

with open('reg.json', 'rt', encoding='utf-8') as gz_file:
    fasta_dict = json.load(gz_file)

with open("ref.W", 'r') as f:
    lines = f.readlines()




aldf=list()
for line in lines:
    linelis=line.split("\t")
    pa=linelis[6]
    ndori=re.findall(r'[><]',pa)
    numbers = re.findall(r'\d+', pa)
    end=pd.DataFrame(numbers)
    end.index=numbers
    end['node']=end.index
    end['fa'] = end['node'].apply(lambda x: fasta_dict.get(x))
    end['chr'] = linelis[3]
    end['ori'] = ndori
    end['length'] = end['fa'].str.len()
    end['end'] = end['length'].cumsum()
    end['start']=end['end']-end['length']+1
    aldf.append(end)

df=(pd.concat(aldf))[['chr','start','end','node','ori']]
df.to_csv("reg.posCHM.txt", sep="\t", index=False, header=False)
