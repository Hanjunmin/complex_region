import argparse

# parser = argparse.ArgumentParser(description="****")

# parser.add_argument("--region", help="***")
# args = parser.parse_args()

# regionuniq=args.region
# regionuniq="chr17:544001-762500"
# regionuniq="chr10:59948501-60002500"
import pandas as pd
import numpy as np
alPath=pd.read_csv("path.txt",sep="\t",header=None)
CHMPath=pd.read_csv("CHM13.path",sep="\t",header=None)
alPath=pd.concat([alPath, CHMPath], axis=0)
alPath['region']=alPath[1].str.split("@|_0_").str[0]
alPath['sam']=alPath[1].str.split("@|_0_").str[1]
alPath[0]=alPath[2].str.count(',')
df_unique = alPath.groupby('region').apply(lambda x: x.drop_duplicates(subset=['sam'])).reset_index(drop=True)
#deleted_rows = alPath.merge(df_unique, how='left', indicator=True).loc[lambda x: x['_merge'] == 'left_only']


alS = pd.read_csv("Sregional.node", sep="\t", header=None, usecols=[1,2])
alS.columns=["S","region"]
alS['S'] = alS['S'].astype(str)


from collections import Counter

for reg in set(df_unique['region']):
    print(reg)
    subpathal=alPath[alPath['region']==reg]
    subdf_unique=df_unique[df_unique['region']==reg]
    subalS=alS[alS['region']==reg]
    seleal=[]
    selectlis={}
    selectlisal={}
    for index,row in subdf_unique.iterrows():
        nodelist = (row[2]).replace('+', '').replace('-', '').split(",")
        nodelist=[i for i in nodelist if i !=""]
        x=nodelist[0:10]
        y=nodelist[-10:]
        seleal.extend(nodelist[0:10])
        seleal.extend(nodelist[-10:])
        selectlis[index]=(x,y)
        selectlisal[index]=x+y
    frequency = Counter()
    for key, (list1, list2) in selectlis.items():
        sorted_tuple = tuple(sorted([tuple(sorted(list1)), tuple(sorted(list2))]))
        frequency[sorted_tuple] += 1
    subalS=subalS[subalS['S'].isin(seleal)]
    for index,value in selectlisal.items():
        df = pd.DataFrame((set(value)), columns=["S"])
        df[index]=1
        subalS = pd.merge(subalS, df, on='S', how='left')
    alSfill = subalS.replace(np.nan, 0)
    row_sums = alSfill.iloc[:, 2:].sum(axis=1)
    alSfill['sum']=row_sums
    alSfill_sorted = alSfill.sort_values(by='sum', ascending=False)
    alSfill_sorted=alSfill_sorted.iloc[:20,]
    Sset=set(alSfill_sorted['S'])
    alist=[]
    for index,row in alSfill_sorted.iterrows():
        node=row['S']
        for i,j in frequency.items():
            if node in i[0]:
                ana=i[1]
            if node in i[1]:
                ana=i[0]
            if len(set(ana) & Sset)!=0:
                Bregion=alSfill_sorted[alSfill_sorted['S'].isin(ana)].iloc[0,]
                alist.append([row['S'],row['sum'],Bregion['S'],Bregion['sum']])
    end=pd.DataFrame(alist)
    if len(end)==0:
        end.to_csv(reg+'.path.dou.filt',index=False,header=False,sep="\t")
        subpathal.to_csv(reg+'.path.dou.bef',index=False,header=False,sep="\t")
        continue
    end['sum']=end[1]+end[3]
    df_sorted = end.apply(lambda row: sorted([row[0], row[2]]), axis=1, result_type='expand')
    df_sorted['sum'] = end['sum']
    grouped = df_sorted.groupby([0,1]).sum().reset_index()
    grouped_sorted = grouped.sort_values(by='sum', ascending=False)
    A=grouped_sorted.iloc[0][0]
    B=grouped_sorted.iloc[0][1]
    print(A)
    print(B)
    selec=[]
    for index,row in subdf_unique.iterrows():
        nodelist = (row[2]).replace('+', '').replace('-', '').split(",")
        nodelist=[i for i in nodelist if i !=""]
        x=nodelist[0:10]
        y=nodelist[-10:]
        if A in x and B in y:
            selec.append(index)
        if A in y and B in x:
            selec.append(index)
    subdf_unique.loc[list(set(selec))].to_csv(reg+'.path.dou.filt',index=False,header=False,sep="\t")
    subpathal.to_csv(reg+'.path.dou.bef',index=False,header=False,sep="\t")

# command='cat '+regionuniq+'.path.dou.end |grep -E "' + str(A)+'+|' + str(A) +'-" |grep -E "' +str(B)+'+|' + str(B) + '-" >'+regionuniq+'.path.dou.end.filt'
# command='cat '+regionuniq+'.path.dou |grep -E "' + str(A)+'+|' + str(A) +'-" |grep -E "' +str(B)+'+|' + str(B) + '-" >'+regionuniq+'.path.dou.filt'

# print(command)
# subprocess.run(command,shell=True)