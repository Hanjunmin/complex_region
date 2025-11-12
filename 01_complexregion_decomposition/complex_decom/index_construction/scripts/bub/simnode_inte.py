import argparse

parser = argparse.ArgumentParser(description='simnode inte file')
parser.add_argument('--reg', help='region')
args = parser.parse_args()


region=args.reg

S_dict={}
P_lis=[]
S_lis=[]


with open(region+'.sim.gfa', 'r') as file:
    for line in file:
        linelis=line.split("\t")
        if linelis[0]=="S":
            S_dict[linelis[1]]=linelis[2]
            S_lis.append([linelis[1],linelis[2]])
        if linelis[0]=="P":
            P_lis.append("|".join(linelis[2].replace("+","").replace("-","").split(",")))


import sys
if len(S_lis)==0:
    simnode=0
    with open('simnode'+region+'.txt', 'w', encoding='utf-8') as f:
        f.write(str(simnode))
    sys.exit()


import pandas as pd
Sdata = pd.DataFrame(S_lis)
Sdata['len']=Sdata[1].str.len()


lensum=sum(Sdata['len'])
Sdata['ratio']=Sdata['len']/lensum
Sdata_sorted = Sdata.sort_values(by='ratio', ascending=False)


results = []
cumulative_sum = 0
for index, row in Sdata_sorted.iterrows():
    cumulative_sum += row["ratio"]
    results.append(row)
    if cumulative_sum > 0.9:
        break

Sdata=pd.DataFrame(results)
# Sdata=Sdata[Sdata['ratio']>0.05]
P_lis=set(P_lis)
savenode=list(Sdata[0])

out=[]
for P in P_lis:
    Pl=P.split("|")
    Pl=[i for i in Pl if i in savenode]
    out.append("|"+("|".join(Pl))+"|")

intelisdict={}
for i in range(len(out)):
    intelisdict[i]=out[i]

intelisdictrev={}
for key,value in intelisdict.items():
    intelisdictrev[key]=('|'.join(str.split(value,"|")[::-1]))



expnode=[]
for m in out:
    elements=str.split(m,"|")
    nodes = [elements[i] +"|"+ elements[i + 1] for i in range(len(elements) - 1)]
    expnode.extend(nodes)
    expnode=list(set(expnode))

unique_elements = set()
for item in expnode:
    parts = sorted(item.split("|"))  # 按字母顺序排序
    unique_elements.add("|".join(parts))  # 重新组合并添加到集合中


savepair=[]
for A in unique_elements:
    L="T"
    Alis=A.split("|")
    vali=set([key for key, value in intelisdict.items() if value.find("|"+A+"|")!= -1]+[key for key, value in intelisdictrev.items() if value.find("|"+A+"|")!= -1])
    subno= intelisdict.keys()- vali
    for sub in subno:
        if intelisdict[sub].find(Alis[0]+"|")!=-1:
            L="S"
            break
    for sub in subno:
        if intelisdict[sub].find(Alis[1]+"|")!=-1:
            L="S"
            break
    if L=="T":
        savepair.append(A)

cannode=[]
for subsave in savepair:
    sublis=subsave.split("|")
    cannode.extend(sublis)

cannode=list(set(cannode))
from collections import Counter
intelis=out
inte=[] ##以每条为基准 至少有20条路径支持，且剩下的一条路径中不存在相关的任何节点
for m in intelis:
    print("#####")
    print(m)
    ben=str.split(m,"|")
    s=0
    e=s+2
    havall=[]
    while e<=len(ben):
        now=ben[s:e]
        letter='|'.join(now)+"|"
        letter_list=set(str.split(letter,"|"))
        letter_list={item for item in letter_list if item}
        if set(now).issubset(set(cannode)):
            vali=set([key for key, value in intelisdict.items() if value.find(letter)!= -1]+[key for key, value in intelisdictrev.items() if value.find(letter)!= -1])
            ###多个segment的必须都要有对应的"pair"
            valinew=[]
            for key in vali:
                value=intelisdict[key]
                valuli=value.split("|")
                list_counter = Counter(valuli)
                appearcount=max({element: list_counter[element] for element in set(letter_list)}.values())
                strcount=value.count(letter)+value.count("|".join([item for item in letter.split("|")[::-1] if item])+"|")
                if appearcount==strcount:
                    valinew.append(key)
            vali=set(valinew)
            subno= intelisdict.keys()- vali
            nocorlis=[]
            for subnohap in subno:
                A=set(str.split(intelisdict[subnohap],"|"))
                if len(A&letter_list)==0:
                    nocorlis.append(subnohap)
            if len(vali)+len(nocorlis)==len(intelisdict.keys()):
                havall.append(letter)
                e=e+1
            else:
                if len(havall)!=0:
                    inte.append(havall[-1])
                s=e-1
                e=s+2
        else:
            s=e-1
            e=s+2
        if len(havall)!=0:
            inte.append(havall[-1])

def filter_inte(df):
    dela=[]
    for sub in df:
        for sub1 in df:
            if sub==sub1:
                continue
            else:
                lar="|"+sub1+"|"
                sm="|"+sub+"|"
                if lar.find(sm)!=-1:
                   dela.append(sub) 
    dfe=[x for x in df if x not in dela]
    return dfe

inte=[s.strip("|") for s in inte]
inte=filter_inte(list(set(inte)))
seen=set()
for item in inte:
    if item not in seen and "|".join(item.split("|")[::-1]) not in seen :
        seen.add(item)

inte=list(seen)


deli=[]
for i in range(0,len(inte)):
    for j in range(i+1,len(inte)):
        if(("|".join(str.split(inte[j],"|")[::-1]))==inte[i]):
            deli.append(inte[j])

inte = [x for x in inte if x not in deli]

Sdata['group']=-1
id=0
for i in inte:
    ilis=i.split("|")
    Sdata.loc[Sdata[0].isin(ilis), 'group'] = id
    id=id+1


for index,row in Sdata.iterrows():
    if row['group']==-1:
        Sdata.loc[index,"group"]=id
        id=id+1


simnode=len(set(Sdata['group']))
with open('simnode'+region+'.txt', 'w', encoding='utf-8') as f:
    f.write(str(simnode))


'''
other_indict_fa={}
other_indict_falis={}
intenum=0
for intenode in list(inte):
    nname="INTEa"+str(intenum)
    intelis = [s.replace(intenode, nname) for s in intelis]
    intelis = [s.replace(("|".join(str.split(intenode,"|")[::-1])), nname) for s in intelis]
    other_indict_falis[nname]=intenode
    intenum=intenum+1

'''