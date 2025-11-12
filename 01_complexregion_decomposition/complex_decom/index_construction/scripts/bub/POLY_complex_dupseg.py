####这个文件处理dup的时候不去split再聚类 这样快很多  我只改了dupiter和dupgene_factlen函数

#没有加anchor的情况

import subprocess
#from Bio import SeqIO
import re
from multiprocessing import Pool
import time
import argparse
parser = argparse.ArgumentParser(description='MERGR')
# parser.add_argument('--p', help='number of process')
# parser.add_argument('--r', help='Reference fasta')
# parser.add_argument('--q', help='query fasta')
parser.add_argument('--sim', help='similarity')
parser.add_argument('--mb', help='mergeblock')
parser.add_argument('--dell', help='delNODElength')
parser.add_argument('--region', help='delNODElength')
args = parser.parse_args()
simarg=float(args.sim)
mblock=int(args.mb)
dellen=int(args.dell)
region=args.region
# simarg=0.9
# mblock=250
# dellen=200
#### 试一下划分窗口然后得到大片段的相似性？？
import numpy as np
from sklearn.cluster import DBSCAN
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

# def get_pa_dictlen(dic):
# allis=[]
# for key,value in padict.items():
#     val=padict[key].replace("+","").replace("-","").split(",")
#     nodefa=list(map(fasta_dict.get, set(val))) 
#     lenn=sum([len(sublist) for sublist in nodefa])
#     allis.append([key,lenn])


def inte_initnode(intelistn,path,nodic):  
    '''
    1.{'INTE0': 'NODE11|NODE12', 'INTE1': 'NODE5|NODE21|NODE22'...}
    2.{'NODE0|+|NODE25|+|NODE1|+|NODE2|+|NODE3|+|NODE4|+|NODE6|+|NODE7|+|NODE8|+|NODE9|+|NODE10|+|NODE11|+|NODE12|+|NODE13|+|NODE23|+|NODE14|+|NODE15|+|NODE16|+|NODE17|+|NODE18|+|NODE19|+|+|NODE20|+|', 'NODE0|+|...
    3.{'NODE0':'114552+,22550220'...}
    '''
    # intelistn=other_indict_falis
    # path=set(rowdictLETTERlist)
    # nodic=padict_letter
    intepath={}
    for name,j in intelistn.items():
        corori=[]
        for m in path: ##循环每个路径寻找
            lis=m.split("|")
            lisnode=[]
            nodeori=[]
            for i in  range(len(lis)-1):
                if lis[i] != "+" and lis[i] != "-":
                    lisnode.append(lis[i])
                    nodeori.append(lis[i+1])
                else:
                    continue
            back="|".join(lisnode)
            pos=back.find(j)
            if pos!=-1:
                pos1=back[:pos].count("|")
                end1=pos1+len(j.split("|"))
                alstr=""
                for loc in list(range(pos1,end1)):
                    corori.append(nodeori[pos1:end1])
                    alstr=alstr+nodic[lisnode[loc]][:-1]+nodeori[loc]+","
                    intepath[name]=alstr
                break
    allpdata1 = pd.DataFrame({'path': list(map(intepath.get, intepath.keys()))})
    allpdata1['h']="P"
    allpdata1['anno']=intepath.keys()
    allpdata1['end']="*"
    return allpdata1

def process(pa):
    palis=pa.split("|")
    pali=[]
    pairlis=[]
    setlis=[]
    for i in range(0,len(palis)-1):
        if palis[i]!="+" and palis[i]!="-":
            cho=palis[i]+palis[i+1]
            pali.append(cho)
            pairlis.append(palis[i]+"|"+palis[i+1])
            setlis.append(palis[i])
    return pali,pairlis,setlis

# def remer(sim_matrix_node):




def inte_p(P,inte):
    allP_dictout={}
    interev={}
    alnode=[]
    for key,value in inte.items():
        interev[key+"_rev"]="|".join(value.split("|")[::-1])
    inte.update(interev)
    for subpkey,subpvalue in P.items():
        can=subpvalue.replace("-","").replace("+","").replace(",","|")
        pah=subpvalue.replace(",","|")
        subtituedict={}
        for intename,intenode in inte.items():
            s=can.find(intenode)
            if s==-1:
                continue
            intenodeid=can[:s].count("|") ###第几个节点开始被inte
            intenodeide=intenodeid+intenode.count("|")
            init="|".join((subpvalue.split(",")[intenodeid:(intenodeide+1)]))
            subtituedict[init[:-1]]=intename
        for key,value in subtituedict.items():
            pah=pah.replace(key, value) 
        allP_dictout[subpkey]=pah.replace("|",",").replace("_rev","") ###这个rev的替换是因为好像重复了。。。之后改一下看看怎么搞
    ###P-----L
    pairsal=[]
    for val in allP_dictout.values():
        val=val.replace("+",",+").replace("-",",-")
        palis=val.split(",")
        pairlis=[]
        for i in range(0,len(palis)-1):
            if palis[i]!="+" and palis[i]!="-":
                cho=palis[i]+palis[i+1]
                pairlis.append(palis[i]+"|"+palis[i+1])     
        pairs = [[pairlis[i],pairlis[i+1]] for i in range(len(pairlis) - 1)]
        pairsal.extend(pairs)
    allL = [[part for item in sublist for part in item.split('|')] for sublist in pairsal]
    return allP_dictout,allL


# def cal_sim():
#     padict_letter





def jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union != 0 else 0

def dict_keys_jaccard_similarity(my_dict):
    keys1 = list(my_dict.keys())
    keys=[]
    for key in keys1:
        if my_dict[key]!="":
            keys.append(key)
    similarity_matrix = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            key1, key2 = keys[i], keys[j]
            val1=my_dict[key1].replace("+","").replace("-","").split(",")
            val2=my_dict[key2].replace("+","").replace("-","").split(",")
            val1nodefa=list(map(fasta_dict.get, set(val1))) 
            val2nodefa=list(map(fasta_dict.get, set(val2))) 
            inter=set(val1)& set(val2)
            internodefa=list(map(fasta_dict.get, inter)) 
            len1=sum([len(sublist) for sublist in val1nodefa])
            len2=sum([len(sublist) for sublist in val2nodefa])
            leninter=sum([len(sublist) for sublist in internodefa])
            similarity = jaccard_similarity(set(val1), set(val2))
            if similarity>0:
                similarity_matrix.append([key1, key2,similarity,leninter,len1,len2])
    return similarity_matrix


# dict_keys_jaccard_similarity(padict)
# out=dict_keys_jaccard_similarity(padict_letter)








import seaborn as sns

# p_data=allpdatanode_seg
# dic=other_indict_falis
# aln=allnode  ##总的allnode

def gene_color(p_data,dic,an,aln):
    p_data['have']=0
    ###group
    if an==0:    ###为0说明颜色代表的是不同node
        p_data['group']=range(len(p_data))
        p_data['have']=1
    else:###颜色代表不同的inte
        for grou in aln:
            if grou in  list(p_data['anno']):
                p_data.loc[p_data['anno'] == grou, 'have'] = 1
                p_data.loc[p_data['anno'] == grou, 'group'] = grou
            else:
                if grou in dic.keys():
                    suli=dic[grou].split("|")
                    p_data.loc[p_data['anno'].isin(suli), 'have'] = 1
                    p_data.loc[p_data['anno'].isin(suli), 'group'] = grou
    p_data=p_data[p_data['have']==1]
    unique_groups = set(p_data['group'])
    colors = sns.color_palette("husl", len(unique_groups))
    hex_colors = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b in colors]
    group_color_dict = dict(zip(unique_groups, hex_colors))
    colordata=p_data.copy()
    colordata['color'] = colordata['group'].map(group_color_dict)
    colorlis=[]
    for i in range(len(colordata)):
        for j in (colordata.iloc[i,]['path'].replace("+","").replace("-","")).split(","):
            colorlis.append([j,colordata.iloc[i,6]])
    # print(pd.DataFrame(colorlis))
    # print(set(pd.DataFrame(colorlis)[1]))
    pd.DataFrame(colorlis).to_csv(region+".color.csv",index=False,header=0,sep=",")
    colordata.to_csv(region+".colorinit.csv",index=False,header=0,sep="\t")

def subgene_factlen(bennode,nownodeall,internode,keysame): ##bench的这个区域的node,某个样本的全部node，有交集的node
    if len(internode)==0:
        return "CONTINUE"
    inter_in_num1_positions = [i for i, x in enumerate(bennode) if x in internode]
    arr_reshaped = np.array(inter_in_num1_positions).reshape(-1, 1)
    dbscan = DBSCAN(eps=20, min_samples=1)  ###如果位置差异大于20就分为另一类
    label=dbscan.fit_predict(arr_reshaped)
    most_common_label = max(set(label.tolist()), key=label.tolist().count)
    res = np.array(inter_in_num1_positions)[label == most_common_label]
    num0=bennode[min(res):(max(res)+1)]
    resfa=list(map(fasta_dict.get, set(num0)))  #end.isnull().values.any()
    lenresfa=sum([len(sublist) for sublist in resfa]) ##query对应ref窗口下真正的长度
    lenresfainit=sum([len(sublist) for sublist in list(map(fasta_dict.get, set(bennode)))]) ##query对应ref窗口下真正的长度
    if len(bennode) == len(internode):
        linkagesim=1
    else:
        linkagesim=(lenresfa)/lenresfainit
    if linkagesim>simarg:
        return "INTE"
    else:
        return "AJUST"

def ajust_region(spli,start,end,rowdict):
    step = (end - start) // spli
    slices = [(i, i + step) for i in range(start, end, step)]
    spliceall=[]
    for sli in slices:
        aduct=[]
        for key12,value12 in rowdict.items():
            num1=bech[sli[0]:sli[1]]
            nowcho=set(value12)
            nextcho=set(num1)
            inter=nowcho & nextcho
            anno=subgene_factlen(num1,value12,inter,key12)
            if anno=="AJUST":
                spliceall.append(1)
                break
    # print(len(spliceall))
    if len(spliceall)/spli>0.5:
        return "POLYCOMPLEX"
    else:
        return "YES"

polycomplex=0




def ajust_regionnew(spli,start,end,rowdict):
    step = (end - start) // spli
    slices = [(i, i + step) for i in range(start, end, step)]
    spliceall=[]
    aduct=[]
    delnum=[]
    for sli in slices:
        for key12,value12 in rowdict.items():
            num1=bech[sli[0]:sli[1]]
            nowcho=set(value12)
            nextcho=set(num1)
            inter=nowcho & nextcho
            anno=subgene_factlen(num1,value12,inter,key12)
            if anno=="AJUST":
                spliceall.append(1)
                aduct.append(sli)
                delnum.extend(num1)
                break
    return aduct,delnum


# bennode=num1
# nownodeall=value12
# internode=inter12
# keysame=key12




def gene_factlen(bennode,nownodeall,internode,keysame): ##bench的这个区域的node,某个样本的全部node，有交集的node
    if len(internode)==0:
        return "CONTINUE",0,0,0,0,keysame
    inter_in_num1_positions = [i for i, x in enumerate(bennode) if x in internode]
    inter_in_num2_positions = [i for i, x in enumerate(nownodeall) if x in internode]
    splitpos = [i for i, x in enumerate(nownodeall) if "_" in x or "dup" in x]
    combined = np.concatenate((inter_in_num2_positions, splitpos))
    sorted_combined = np.sort(combined)
    combinend = np.where(np.isin(sorted_combined, splitpos), '*', sorted_combined)
    splilis="_".join(combinend).split("*")
    maxelen=[]
    for element in splilis:
        temp=element.split("_")
        temp=[int(float(i)) for i in temp if i != '']
        if len(temp)>len(maxelen):
            maxelen=temp
    inter_in_num2_positions=maxelen
    nseg=[nownodeall[i] for i in maxelen]
    internode=set(nseg) & internode
    inter_in_num1_positions = [i for i, x in enumerate(bennode) if x in internode]
    arr_reshaped = np.array(inter_in_num1_positions).reshape(-1, 1)
    arr_reshaped2 = np.array(inter_in_num2_positions).reshape(-1, 1)
    dbscan = DBSCAN(eps=40, min_samples=1)  ###如果位置差异大于20就分为另一类
    label=dbscan.fit_predict(arr_reshaped)
    label2=dbscan.fit_predict(arr_reshaped2)
    most_common_label = max(set(label.tolist()), key=label.tolist().count)
    most_common_label2 = max(set(label2.tolist()), key=label2.tolist().count)
    res = np.array(inter_in_num1_positions)[label == most_common_label]
    res2 = np.array(inter_in_num2_positions)[label2 == most_common_label2]
    num0=nownodeall[min(res2):(max(res2)+1)]
    num0=[item for item in num0 if item not in ['+', '-','']]
    # print(set(num0))
    num0fa=list(map(fasta_dict.get, set(num0)))  #end.isnull().values.any()
    num0fa = [item for item in num0fa if item is not None]
    # print(num0fa)
    lenquef=sum([len(sublist) for sublist in num0fa]) ##query对应ref窗口下真正的长度
    if len(bennode) == len(internode):
        linkagesim=1
        queal=nownodeall[min(res2):(max(res2)+1)]
        queal=[item for item in queal if item not in ['+', '-','']]
        quealfa=list(map(fasta_dict.get, set(queal)))  #end.isnull().values.any()
        quealfa = [item for item in quealfa if item is not None]
        lenque=sum([len(sublist) for sublist in quealfa]) ##query对应ref窗口下真正的长度
        linkagequesim=lenque/lenquef  ###query对应的inter的最大和最小的node对应的query的长度
        simqueref=linkagequesim
    else:
        benal=bennode[min(res):(max(res)+1)]
        queal=nownodeall[min(res2):(max(res2)+1)]
        queal=[item for item in queal if item not in ['+', '-','']]
        benalfa=list(map(fasta_dict.get, set(benal)))  #end.isnull().values.any()
        benalfa = [item for item in benalfa if item is not None]
        quealfa=list(map(fasta_dict.get, set(queal)))  #end.isnull().values.any()
        quealfa = [item for item in quealfa if item is not None]
        lenben=sum([len(sublist) for sublist in benalfa]) ##query对应ref窗口下真正的长度
        lenque=sum([len(sublist) for sublist in quealfa]) ##query对应ref窗口下真正的长度
        x=list(map(fasta_dict.get, set(bennode)))
        x=[item for item in x if item is not None]
        lenbeninit=sum([len(sublist) for sublist in x]) ##query对应ref窗口下真正的长度
        linkagesim=max(((max(res)-min(res))+1)/len(bennode),lenben/lenbeninit)
        linkagequesim=lenque/lenquef  ###query对应的inter的最大和最小的node对应的query的长度
        simqueref=min(lenque,lenbeninit)/max(lenque,lenbeninit)
    # print(linkagesim)
    # print(linkagequesim)
    # posal=[pos for pos, val in enumerate(nownodeall) if val in internode]
    # arr_reshaped = np.array(posal).reshape(-1, 1)
    # dbscan = DBSCAN(eps=(len(internode)), min_samples=1)  ###from sklearn.cluster import DBSCAN
    # label=dbscan.fit_predict(arr_reshaped)
    # for lab in set(label):
    #     selected_rows = [row for row, l in zip(posal, label) if l == lab]
    selected_rows = res2
    posma=max(selected_rows)  
    posmin=min(selected_rows)
    truenode=nownodeall[posmin:(posma+1)] ###真正的路径
    truenode=[item for item in truenode if item not in ['+', '-','']]
    truenodefa=list(map(fasta_dict.get, truenode)) 
    if None in truenodefa:
        return "CONTINUE",0,0,0,0,keysame
    lena=sum([len(sublist) for sublist in truenodefa]) ##query对应ref窗口下真正的长度
    nextchofa=list(map(fasta_dict.get, set(bennode)))  #end.isnull().values.any()
    interfa=list(map(fasta_dict.get, internode))  ####有交集的路径
    lenb=sum([len(sublist) for sublist in nextchofa]) ##ref真正的长度
    leninter=sum([len(sublist) for sublist in interfa]) ##有交集真正的长度
    fro=rowdictalpath[keysame].find(nownodeall[posmin])  
    to=rowdictalpath[keysame].find(nownodeall[posma])+len(nownodeall[posma])+1
    # print(lena)
    # print(lenb)
    # print(leninter)
    # print(min(lena,lenb)/max(lena,lenb))
    # print(leninter/lenb)
    # print(linkagesim)
    # if min(lena,lenb)/max(lena,lenb)>simarg and leninter/lenb>simarg:
    if linkagesim>simarg and simqueref>simarg:
        return "INTE",fro,to,posmin,posma,keysame
    else:
        return "AJUST",fro,to,posmin,posma,keysame 






def iter(padict,chunk_size,rowdict,bech,bechname,repid,delnode,polycomplex): ###winow,全部的path，某样本大于500的path（bench）
    start=0
    while start<len(bech):
        if (len(bech)<chunk_size*1.5) and len(bech)>chunk_size: ###说明chunk要比bech小一点点 如果按照chunk来会损失一些node或者报错 
            chunk_size=len(bech)
        end=start+chunk_size
        na=bechname+"_"+str(start)+"_"+str(end)+"@"+str(repid)
        print(start,end)
        annolis=[]
        num1=bech[start:end]
        numbef=num1
        ###将该区域分成spli份 判断出现AJUST的次数  如果出现次数是spli的0.5以上 则规定这段区域为一个多态的区域 将node后面加_并且count+1为多态的系数
        na1,na2=ajust_regionnew(5,start,end,rowdict)
        print(na1)
        print(na2)
        if len(na1)>2: ###感觉50个segment比较合适
            replace_map = {str(item): f'_{item}' for item in na2}
            rowdict = {
                k: [replace_map.get(item, item) for item in v]
                for k, v in rowdict.items()
            }
            delnode.extend(na2)
            print("delanno:")
            print(na2)
            delnode=list(set(delnode))
            start=end
            polycomplex=polycomplex+1
            print("!!!!!!!!!!!!!!!!!!!!!")
            print(polycomplex)
            continue
        dele="false"
        for key12,value12 in rowdict.items():
            nowcho12=set(value12)
            nextcho12=set(num1)
            inter12=nowcho12 & nextcho12
            [anno,fro,to,posmin,posma,keya]=gene_factlen(num1,value12,inter12,key12)
            if anno=="INTE":
                annolis.append(["INTE",start,end,keya])
            if anno=="AJUST":
                positions = [num1.index(x) for x in inter12 if x in num1] ##交集的node在原始的序列中的位置 根据这个位置进行迭代直到差异很大为止
                ne=0
                rev=sorted(positions)
                s=start+min(rev)
                eid=len(rev)-1
                while eid>0:
                    e=start+rev[eid]+1
                    num2=bech[s:e]
                    nextcho1=set(num2)
                    inter1=nowcho12 & nextcho1
                    [annoa,fro,to,posmin,posma,keya]=gene_factlen(num2,value12,inter1,key12)
                    if annoa=="INTE":
                        ne=e
                        break
                    else:
                        eid=eid-1
                if ne!=0:
                    annolis.append(["AJUST",start+min(sorted(positions)),ne,keya])
                if ne==0 and len(inter12)/len(num1)>0.5:
                    dele="true"
                    break
        if dele=="true":
            start=end
            continue
        if len(annolis)==0: ###说明这段区域在不同样本中变化很大 左右扩展边界也达不到相似性大于0.9
            print("Oh no")
            print(num1)
            replace_map = {str(item): f'_{item}' for item in num1}
            rowdict = {
                k: [replace_map.get(item, item) for item in v]
                for k, v in rowdict.items()
            }
            delnode.extend(num1)
            print("delanno:")
            print(num1)
            delnode=list(set(delnode))
            start=end
            continue
        print("yes")
        allhap=pd.DataFrame(annolis)
        start=max(allhap[1])
        end=min(allhap[2])
        na=bechname+"_"+str(start)+"_"+str(end)+"@"+str(repid)
        print("yes1")
        num1=bech[start:end]
        lennum1=len("".join(list(map(fasta_dict.get,num1))))
        if lennum1<dellen:
            replace_map = {str(item): f'_{item}' for item in num1}
            rowdict = {
                k: [replace_map.get(item, item) for item in v]
                for k, v in rowdict.items()
            }
            delnode.extend(num1)
            print("delanno:")
            print(num1)
            delnode=list(set(delnode))
            start=end
        else:
            for keyx in set(list(allhap[3])):
                valuex=rowdict[keyx]
                nowcho12=set(valuex)
                inter12=nowcho12 & set(numbef)
                nowcho=set(valuex)
                nextcho=set(num1)
                inter=nowcho & nextcho
                [anno,fro,to,posmin,posma,keya]=gene_factlen(num1,valuex,inter,keyx)
                valuex[posmin:(posma+1)]=[na]
                rowdict[keyx]=valuex
                # print("bbbbbbbbbbbbbbbbbbbb")
                # if keyx=='C023-CHA-S03_2-Pat_0_C023-CHA-S03_2_Pat_chr22_14654170_46858844':
                #     print(rowdict[keyx])
                # print("aaaaaaaaaaaaaaaaaaa")
                padict[na]=rowdictalpath[keyx][fro:to]
                # print(na)
                # print(keyx)
            start=end
            # print(rowdict['C023-CHA-S03-Pat_0_C023-CHA-S03_Pat_chr22_18409025_48958486'])
            print(padict)
    return(rowdict,padict,delnode,polycomplex)

# bennode=num1
# nownodeall=value12
# internode=inter12
# keysame=key12

def dupgene_factlen(bennode,nownodeall,internode,keysame,start,end): 
    if len(internode)==0:
        return [["CONTINUE",0,0,0,0,0,0,0,0,0]]
    x=list(map(fasta_dict.get, set(bennode)))
    x=[item for item in x if item is not None]
    lenbeninit=sum([len(sublist) for sublist in x]) ##query对应ref窗口下真正的长度
    inter_in_num1_positions = [i for i, x in enumerate(bennode) if x in internode]
    inter_in_num2_positions = [i for i, x in enumerate(nownodeall) if x in internode]
    splitpos = [i for i, x in enumerate(nownodeall) if "_" in x]
    combined = np.concatenate((inter_in_num2_positions, splitpos))
    sorted_combined = np.sort(combined)
    combinend = np.where(np.isin(sorted_combined, splitpos), '*', sorted_combined)
    splilis="_".join(combinend).split("*")
    maxelenn=[]
    # print(splilis)
    for element in splilis:
        temp=element.split("_")
        temp=[int(float(i)) for i in temp if i != '']
        if len(temp)>len(maxelenn) and len(temp)>2:
            maxelenn=temp
    if len(maxelenn)==0:
        return [["CONTINUE",0,0,0,0,0,0,0,0,0]]
    # print(maxelenn)
    inter_in_num2_positions=[]
    label2=[]
    ext=1
    for element in splilis:
        temp=element.split("_")
        temp=[int(float(i)) for i in temp if i != '']
        if len(temp)/len(maxelenn)>0.5:
            inter_in_num2_positions.extend(temp)
            arr_reshaped2 = np.array(temp).reshape(-1, 1)
            dbscan = DBSCAN(eps=100, min_samples=1)  ###如果位置差异大于20就分为另一类
            labeln=dbscan.fit_predict(arr_reshaped2)
            labeln=labeln+ext
            ext=max(labeln)+1
            label2.extend(labeln)
    # print(inter_in_num2_positions)
    # print(label2)
    nseg=[nownodeall[i] for i in inter_in_num2_positions]
    internode=set(nseg) & internode
    outlis=[]
    for lab in set(label2):
        res2 = np.array(inter_in_num2_positions)[np.array(label2) == lab]
        num0=nownodeall[min(res2):(max(res2)+1)]
        num0=[item for item in num0 if item not in ['+', '-','']]
        num0fa=list(map(fasta_dict.get, set(num0)))  #end.isnull().values.any()
        internodesub=set(num0) & internode
        inter_in_num1_positions = [i for i, x in enumerate(bennode) if x in internodesub]
        arr_reshaped = np.array(inter_in_num1_positions).reshape(-1, 1)
        label=dbscan.fit_predict(arr_reshaped)
        most_common_label = max(set(label.tolist()), key=label.tolist().count)
        res = np.array(inter_in_num1_positions)[label == most_common_label]
        if len(bennode) == len(internodesub):
            linkagesim=1
            queal=nownodeall[min(res2):(max(res2)+1)]
            queal=[item for item in queal if item not in ['+', '-','']]
            quealfa=list(map(fasta_dict.get, set(queal)))  #end.isnull().values.any()
            quealfa = [item for item in quealfa if item is not None]
            lenque=sum([len(sublist) for sublist in quealfa]) ##query对应ref窗口下真正的长度
            linkagequesim=min(lenque,lenbeninit)/max(lenque,lenbeninit)    ###query对应的inter的最大和最小的node对应的query的长度
            simqueref=linkagequesim
        else:
            benal=bennode[min(res):(max(res)+1)]
            queal=nownodeall[min(res2):(max(res2)+1)]
            queal=[item for item in queal if item not in ['+', '-','']]
            benalfa=list(map(fasta_dict.get, set(benal)))  #end.isnull().values.any()
            benalfa = [item for item in benalfa if item is not None]
            quealfa=list(map(fasta_dict.get, set(queal)))  #end.isnull().values.any()
            quealfa = [item for item in quealfa if item is not None]
            lenben=sum([len(sublist) for sublist in benalfa]) ##query对应ref窗口下真正的长度
            lenque=sum([len(sublist) for sublist in quealfa]) ##query对应ref窗口下真正的长度
            linkagesim=max(((max(res)-min(res))+1)/len(bennode),lenben/lenbeninit)
            simqueref=min(lenque,lenbeninit)/max(lenque,lenbeninit) 
        selected_rows = res2
        posma=max(selected_rows)  
        posmin=min(selected_rows)
        truenode=nownodeall[posmin:(posma+1)] ###真正的路径
        truenode=[item for item in truenode if item not in ['+', '-','']]
        truenodefa=list(map(fasta_dict.get, truenode)) 
            # if None in truenodefa:
            #     return "CONTINUE",0,0,0,0
        lena=sum([len(sublist) for sublist in truenodefa]) ##query对应ref窗口下真正的长度
        nextchofa=list(map(fasta_dict.get, set(bennode)))  #end.isnull().values.any()
        interfa=list(map(fasta_dict.get, internodesub))  ####有交集的路径
        lenb=sum([len(sublist) for sublist in nextchofa]) ##ref真正的长度
        leninter=sum([len(sublist) for sublist in interfa]) ##有交集真正的长度
        fro=rowdictalpath[keysame].find(nownodeall[posmin])  
        to=rowdictalpath[keysame].find(nownodeall[posma])+len(nownodeall[posma])+1
        # print(simqueref)
        # print(linkagesim)
        if linkagesim>0.7 and simqueref>0.7:
            outlis.append(["INTE",start,end,fro,to,posmin,posma,lab,keysame,len(res2)])
        else:
            outlis.append(["AJUST",start,end,fro,to,posmin,posma,lab,keysame,len(res2)])
    return outlis


def dupiter(padict,chunk_size,rowdict,bech,bechname,repid,delnode,polycomplex): ###winow,全部的path，某样本大于500的path（bench）
    start=0
    while start<len(bech):
        if (len(bech)<chunk_size*1.5) and len(bech)>chunk_size: ###说明chunk要比bech小一点点 如果按照chunk来会损失一些node或者报错 
            chunk_size=len(bech)
        end=start+chunk_size
        start=int(start)
        end=int(end)
        na=bechname+"_"+str(start)+"_"+str(end)+"@"+str(repid)
        print(start,end)
        annolis=[]
        intesave=[]
        num1=bech[start:end]
        numbef=num1
        ###将该区域分成spli份 判断出现AJUST的次数  如果出现次数是spli的0.5以上 则规定这段区域为一个多态的区域 将node后面加_并且count+1为多态的系数
        na1,na2=ajust_regionnew(5,start,end,rowdict)
        print(na1)
        print(na2)
        if len(na1)>2: ###感觉50个segment比较合适
            replace_map = {str(item): f'_{item}' for item in na2}
            rowdict = {
                k: [replace_map.get(item, item) for item in v]
                for k, v in rowdict.items()
            }
            delnode.extend(na2)
            print("delanno:")
            print(na2)
            delnode=list(set(delnode))
            start=end
            polycomplex=polycomplex+1
            print("!!!!!!!!!!!!!!!!!!!!!")
            print(polycomplex)
            continue
        ajustnum=0
        dele="false"
        for key12,value12 in rowdict.items():
            nowcho12=set(value12)
            nextcho12=set(num1)
            inter12=nowcho12 & nextcho12
            eva=dupgene_factlen(num1,value12,inter12,key12,start,end)
            # [anno,fro,to,posmin,posma]=gene_factlen(num1,value12,inter12,key12)
            if "AJUST" in set(pd.DataFrame(eva)[0]):
                ajustnum=ajustnum+1
                positions = [num1.index(x) for x in inter12 if x in num1] ##交集的node在原始的序列中的位置 根据这个位置进行迭代直到差异很大为止
                ne=0
                rev=sorted(positions)
                s=start+min(rev)
                eid=len(rev)-1
                while eid>0:
                    e=start+rev[eid]+1
                    num2=bech[s:e]
                    nextcho1=set(num2)
                    inter1=nowcho12 & nextcho1
                    evan=dupgene_factlen(num2,value12,inter1,key12,start+min(sorted(positions)),e)
                    if "AJUST" in set(pd.DataFrame(evan)[0]):
                        eid=eid-1
                    if ("CONTINUE" in set(pd.DataFrame(evan)[0])):
                        break
                    if set(pd.DataFrame(evan)[0])=={"INTE"}:
                        ne=e
                        break   
                if ne!=0:
                    annolis.extend(evan)
                if ne==0 and len(inter12)/len(num1)>0.5:
                    dele="true"
                    break
            else:
                annolis.extend(eva)   
        if dele=="true":
            start=end
            continue
        print(annolis)
        print(ajustnum)
        if len(annolis)==0: ###说明这段区域在不同样本中变化很大 左右扩展边界也达不到相似性大于0.9
            print("Oh no")
            print(num1)
            replace_map = {str(item): f'_{item}' for item in num1}
            rowdict = {
                k: [replace_map.get(item, item) for item in v]
                for k, v in rowdict.items()
            }
            delnode.extend(num1)
            print("delanno:")
            print(num1)
            delnode=list(set(delnode))
            start=end
            continue
        print("yes")
        allhap=pd.DataFrame(annolis)
        allhap=allhap[allhap[0]!="CONTINUE"]
        if len(allhap)==0:
            start=end
            continue
        start=max(allhap[1])
        end=min(allhap[2])
        na=bechname+"_"+str(start)+"_"+str(end)+"@"+str(repid)
        print("yes1")
        num1=bech[int(start):int(end)]
        lennum1=len("".join(list(map(fasta_dict.get,num1))))
        if lennum1<dellen:
            replace_map = {str(item): f'_{item}' for item in num1}
            rowdict = {
                k: [replace_map.get(item, item) for item in v]
                for k, v in rowdict.items()
            }
            delnode.extend(num1)
            print("delanno:")
            print(num1)
            delnode=list(set(delnode))
            start=end
        else:
            if ajustnum==0:
                grouped = allhap.groupby(8)
                result = grouped.agg({3: 'min', 4: 'min'}).reset_index()
                result.columns = [8, 'min_5', 'min_6']
                result = allhap.merge(result, on=8, how='left')
                allhap = result.sort_values(by=5, ascending=False)
                for ind,rowx in allhap.iterrows():
                    keyx=rowx[8]
                    valuex=rowdict[keyx]
                    posmin=rowx[5]
                    posma=rowx[6]
                    valuex[posmin:(posma+1)]=[na]
                    rowdict[keyx]=valuex
                    padict[na]=rowdictalpath[keyx][rowx['min_5']:rowx['min_6']]
                start=end
                continue
            for keyx in set(list(allhap[8])):
                valuex=rowdict[keyx]
                nowcho12=set(valuex)
                inter12=nowcho12 & set(numbef)
                nowcho=set(valuex)
                nextcho=set(num1)
                inter=nowcho & nextcho
                # in1df,in1out=splitpath(valuex)
                evanq=dupgene_factlen(num1,valuex,inter,keyx,0,0)
                evanqdf = pd.DataFrame(evanq)
                evanqdf = evanqdf.sort_values(by=5, ascending=False)
                for ind,rowx in evanqdf.iterrows():
                    posmin=rowx[5]
                    posma=rowx[6]
                    valuex[posmin:(posma+1)]=[na]
                rowdict[keyx]=valuex
                padict[na]=rowdictalpath[keyx][min(evanqdf[3]):min(evanqdf[4])]
            start=end
            # print(rowdict['C023-CHA-S03-Pat_0_C023-CHA-S03_Pat_chr22_18409025_48958486'])
            print(padict)
    return(rowdict,padict,delnode,polycomplex)


def split_fasta(fasta_dict,splitbase):
    ####1.split fasta
    # splitbase=10 ###10bp as the split unit 
    fasta_dictsplit={}
    for index,value in fasta_dict.items():
        if len(value)>10:
            splilen=math.floor(len(value) / splitbase)
            fasta_dictsplit[index]=[str(index) for i in range(0, splilen)]
            # fasta_dictsplit[index]=[str(index) +"#"+ str(i) for i in range(0, splilen)]
        else:
            # fasta_dictsplit[index]=[str(index) +"#0"]
            fasta_dictsplit[index]=[str(index)]
    return fasta_dictsplit


def splitpath(in1):
    ###生成每条path对应的segement和对应的split后的多个顶点的位置 in1out in1df
    in1df=[]
    start=0
    in1out = []
    for item in in1:
        if item in keys_set:  # 快速查找
            in1out.extend(fasta_dictsplit[item])
            end=start+len(fasta_dictsplit[item])    
        else:
            in1out.append(item)
            end=start+1
        in1df.append([item,start,end])
        start=end
    in1df=pd.DataFrame(in1df)
    in1df['clus']=-1
    in1df[2]=in1df[2]-1
    return [in1df,in1out]



def subgene_factlen(bennode,nownodeall,internode,keysame): ##bench的这个区域的node,某个样本的全部node，有交集的node
    if len(internode)==0:
        return "CONTINUE"
    inter_in_num1_positions = [i for i, x in enumerate(bennode) if x in internode]
    arr_reshaped = np.array(inter_in_num1_positions).reshape(-1, 1)
    dbscan = DBSCAN(eps=20, min_samples=1)  ###如果位置差异大于20就分为另一类
    label=dbscan.fit_predict(arr_reshaped)
    most_common_label = max(set(label.tolist()), key=label.tolist().count)
    res = np.array(inter_in_num1_positions)[label == most_common_label]
    num0=bennode[min(res):(max(res)+1)]
    resfa=list(map(fasta_dict.get, set(num0)))  #end.isnull().values.any()
    lenresfa=sum([len(sublist) for sublist in resfa]) ##query对应ref窗口下真正的长度
    lenresfainit=sum([len(sublist) for sublist in list(map(fasta_dict.get, set(bennode)))]) ##query对应ref窗口下真正的长度
    if len(bennode) == len(internode):
        linkagesim=1
    else:
        linkagesim=(lenresfa)/lenresfainit
    if linkagesim>simarg:
        return "INTE"
    else:
        return "AJUST"

def ajust_region(spli,start,end,rowdict):
    step = (end - start) // spli
    slices = [(i, i + step) for i in range(start, end, step)]
    spliceall=[]
    for sli in slices:
        aduct=[]
        for key12,value12 in rowdict.items():
            num1=bech[sli[0]:sli[1]]
            nowcho=set(value12)
            nextcho=set(num1)
            inter=nowcho & nextcho
            anno=subgene_factlen(num1,value12,inter,key12)
            if anno=="AJUST":
                spliceall.append(1)
                break
    # print(len(spliceall))
    if len(spliceall)/spli>0.5:
        return "POLYCOMPLEX"
    else:
        return "YES"

import subprocess
#from Bio import SeqIO
import re
from multiprocessing import Pool
import time
import argparse
# import numpy as np
from sklearn.cluster import DBSCAN

import pandas as pd 
import re
import json
from Bio.Seq import Seq
import lmdb
import numpy as np
import math
sdpa=pd.read_csv("temp."+region+".S.start",sep="\t",header=None)
sdpa[0]=sdpa[0].astype(str)
fasta_dict=dict(zip(sdpa.iloc[:, 0], sdpa.iloc[:, 1]))
splitbase=10
fasta_dictsplit=split_fasta(fasta_dict,splitbase)



# with open('/home/jmhan/pan-SD/CHM13-APGp1-HPRCp1-HGSVCp3_MC.S.json', 'rt', encoding='utf-8') as gz_file:
#     fasta_dict = json.load(gz_file)
keys_set = set(fasta_dictsplit.keys())

allpath = pd.read_csv("temp."+region+".allpath",sep="\t",header=None)
rowlist=[]
for index,row in allpath.iterrows():
    rowlist.append(row[1])

rowdictalpath={}
for index,row in allpath.iterrows():
    rowdictalpath[row[1]]=row[2]



from collections import Counter
dupelement=[]
rowdict={}
for index,row in allpath.iterrows():
    in1=str.split(row[2].replace("-",",-").replace("+",",+"),",")
    rowdict[row[1]]=in1  # rowdict[row[1]]=re.findall(r'\d+', row[2])
    counted_elements = Counter(in1)
    for element, count in counted_elements.items():
        if count > 1:
            # print(index)
            dupelement.append(element)






import numpy as np
dupelement=set(dupelement)
dupelement=[i for i in dupelement if i not in ["+","-"]]
dupelemental=[]
###把dup的区域注释一下:
rowdict={}
testlis=[]
apppdata=[]
for index,row in allpath.iterrows():
    print(index)
    in1=str.split(row[2].replace("-",",-").replace("+",",+"),",")
    in1df,in1out=splitpath(in1)
    if len(dupelement)!=0:
        positions = np.where(np.isin(in1out, list(dupelement)))[0].tolist()
        print(len(positions))
        if len(positions)==0:
            rowdict[row[1]]=in1 
            continue
        arr_reshaped = np.array(positions).reshape(-1, 1)
        dbscan = DBSCAN(eps=100, min_samples=1)  ###如果位置差异大于20就分为另一类
        label=dbscan.fit_predict(arr_reshaped)
        for lab in set(label):
            res = np.array(positions)[label == lab]
            in1df.loc[(in1df[1]<=max(res)) & (in1df[2]>=min(res)),'clus']=lab
            xt=in1df.loc[(in1df[1]<=max(res)) & (in1df[2]>=min(res))]
            res=[xt.index[0],xt.index[-1]]
            pa="".join(in1[min(res):(max(res)+1)]).replace("+","+,").replace("-","-,")+"+"
            apppdata.append(["P",str(index)+"@"+str(lab),pa,"*"])
            dupelemental.extend(in1[min(res):(max(res)+1)])
            in1[min(res):(max(res)+1)]=[item + "dup" for item in in1[min(res):(max(res)+1)]]
    dup_count = sum(1 for item in in1 if 'dup' in item)
    rowdict[row[1]]=in1 
    testlis.append([row[1],dup_count/len(in1)])
    dupelemental=list(set(dupelemental))

dupelemental=set(dupelemental)


# pd.DataFrame(apppdata).to_csv("P.data",index=False,header=False,sep="\t")

######non-dup
polycomplex=0
padict={}
# for anchor in rangenode:
#     padict["anchor_"+anchor]=anchor+"+"
delnode=[]
padict={}
repid=0
for subchunk in range(int(mblock),0,-5):
    print(subchunk)
    for key1,value1 in rowdict.items():
        value1=rowdict[key1]
        modified_list = ["*" if "_" in item else item for item in value1]  ###观察连续的元素有没有大于500
        modified_list = ["*" if "dup" in item else item for item in modified_list]
        moti=[x for x in "|".join(modified_list).replace("*|", "*").replace("|*", "*").split("*") if x] 
        reppro= [item for item in moti if item.count("|") > 2*subchunk]
        bechname=key1
        # print(len(reppro))
        if len(reppro)>0: ##说明有segment大于500个的path
            dellisn=[]
            savlisn=[]
            for rep in reppro:
                bech=rep.split("|")
                bech=[item for item in bech if item not in ['+', '-','']]
                if len(set(bech) & dupelemental)/(len(set(bech)))>0.5:
                    replace_map = {str(item): f'{item}dup' for item in bech}
                    rowdict[key1]=[replace_map.get(item, item) for item in value1]
                    continue
                lenbech=len("".join(list(map(fasta_dict.get, bech))))
                print("#########################")
                # print(lenbech)
                print("#########################")
                if lenbech<dellen:
                    dellisn.extend(bech)
                else:
                    savlisn.append(rep)
            replace_map = {str(item): f'_{item}' for item in dellisn}
            rowdict = {
                k: [replace_map.get(item, item) for item in v]
                for k, v in rowdict.items()
            }
            delnode.extend(dellisn)
            delnode=list(set(delnode))
            for rep in savlisn:
                bech=rep.split("|")
                bech=[item for item in bech if item not in ['+', '-','']]
                print("a")
                [rowdict,padict,delnode,polycomplex]=iter(padict,subchunk,rowdict,bech,bechname,repid,delnode,polycomplex)
                print("b")
                # print(rowdict['C023-CHA-S03_2-Pat_0_C023-CHA-S03_2_Pat_chr22_14654170_46858844'])
                repid=repid+1



aaaa=rowdict.copy()
for key1,value1 in rowdict.items():
    value1=rowdict[key1]
    rowdict[key1]=[item.replace('dup', '') if "dup" in item else item for item in value1]


# replace_map = { f'_{item}':str(item) for item in sdpa[0]}
# rowdict = {
#     k: [replace_map.get(item, item) for item in v]
#     for k, v in rowdict.items()
# }
# polycomplex=0
######dup
for subchunk in range(int(mblock),0,-5):
    print(subchunk)
    for key1,value1 in rowdict.items():
        value1=rowdict[key1]
        modified_list = ["*" if "_" in item else item for item in value1]  ###观察连续的元素有没有大于500
        moti=[x for x in "|".join(modified_list).replace("*|", "*").replace("|*", "*").split("*") if x] 
        reppro= [item for item in moti if item.count("|") > 2*subchunk]
        bechname=key1
        # print(len(reppro))
        if len(reppro)>0: ##说明有segment大于500个的path
            dellisn=[]
            savlisn=[]
            for rep in reppro:
                bech=rep.split("|")
                bech=[item for item in bech if item not in ['+', '-','']]
                lenbech=len("".join(list(map(fasta_dict.get, bech))))
                print("#########################")
                # print(lenbech)
                print("#########################")
                if lenbech<dellen:
                    dellisn.extend(bech)
                else:
                    savlisn.append(rep)
            replace_map = {str(item): f'_{item}' for item in dellisn}
            rowdict = {
                k: [replace_map.get(item, item) for item in v]
                for k, v in rowdict.items()
            }
            delnode.extend(dellisn)
            delnode=list(set(delnode))
            for rep in savlisn:
                print("a")
                [rowdict,padict,delnode,polycomplex]=dupiter(padict,subchunk,rowdict,bech,bechname,repid,delnode,polycomplex)
                print("b")
                # print(rowdict['C023-CHA-S03_2-Pat_0_C023-CHA-S03_2_Pat_chr22_14654170_46858844'])
                repid=repid+1

# merrowdict={}
# for index,value in rowdict.items():
#     matches = re.search("Mat|Pat|hap1|hap2|CHM13|CN1|Ch38", index)
#     name=index[:matches.end()]
#     if name in merrowdict:
#         merrowdict[name].extend(value)
#     else:
#         merrowdict[name]=value






# interal=set(merrowdict[list(merrowdict.keys())[0]])
# for index,value in  merrowdict.items():
#     interal= interal&set(value)



import re

pattern = re.compile(r'^_\d+')
rowdict = {
    k: [item for item in v if not pattern.match(item)]
    for k, v in rowdict.items()
}

# print(new_rowdict)
with open('polycomplex'+region+'.txt', 'w', encoding='utf-8') as f:
    f.write(str(polycomplex))

import pickle

data_to_save = (rowdict,padict)


with open('.rowdict'+region+'.pkl', 'wb') as f:
    pickle.dump(data_to_save, f)


import string
keys = list(padict.keys()) 
letter_dict = {key: letter for key, letter in zip(keys, [f"NODE{i}" for i in range(len(keys))])}
padict_letter={}
for key,value in padict.items():
    padict_letter[letter_dict[key]]=padict[key]


with open("out_INTE_node1"+region+".json", "w", encoding="utf-8") as f:
    json.dump(letter_dict, f, ensure_ascii=False, indent=4)


##padict:key-整合的node名字+对应的路径
##rowdictLETTER:
# rowdictLETTER={}
# for key,value in rowdict.items():
#     rowdictLETTER[key]=rowdict['va']


###这里前面加一步 把node比较多的区域和已有的path进行比较
change_dict_delnode={}
for key,value in rowdict.items():
    numbers = [item for item in value if str(item).replace('.', '', 1).isdigit()]
    if len( numbers)==0:
        delnode=[]
    else:
        Snum = pd.DataFrame({'fa': list(map(fasta_dict.get, list(numbers)))})
        Snum['node']=list(numbers)
        Snum['len']=Snum['fa'].str.len()
        delnode=list(Snum[Snum['len']<1000]['node'])
    change_dict_delnode[key] = [item for item in value if item not in delnode]



rowdictLETTER={}
indic={}
rowdictLETTERlist=[]
allnode=[]
for key,value in change_dict_delnode.items():
    rowdictLETTER[key]=[letter_dict.get(key1, key1) for key1 in value]
    rowdictLETTERlist.append('|'.join([letter_dict.get(key1, key1) for key1 in value]))
    indic[key]='|'.join([letter_dict.get(key1, key1) for key1 in value])
    allnode.extend([letter_dict.get(key1, key1) for key1 in value])

allnode=list(set(allnode))
allnode=[item for item in allnode if item not in ['+', '-','']]

allL=[]
allP_dict={}
uniqpath=[]
allnodes_dict={}
id=0
for pa in list(set(rowdictLETTERlist)):
    [a,b,nodes]= process(pa)
    pairs = [[b[i],b[i+1]] for i in range(len(b) - 1)]
    split_list = [[part for item in sublist for part in item.split('|')] for sublist in pairs]
    allL.extend(split_list)
    allP_dict[id]=",".join(a)
    allnodes_dict[id]="|".join(nodes)
    uniqpath.append("|".join(nodes))
    id=id+1

# allL=pd.DataFrame(allL).drop_duplicates()

# allP_dict={}
# id=0
# for pa in rowdictLETTERlist:
#     allnode.extend(str.split(pa,"|"))
#     allP_dict[id]=pa.replace("|","+,")
#     id=id+1

allnode=set(allnode)

import sys
if not padict_letter:
    print("empty node dict")
    sys.exit()

allL=pd.DataFrame(allL).drop_duplicates()
allpdata = pd.DataFrame({'path': list(map(allP_dict.get, allP_dict.keys()))})
allpdata['h']="P"
allpdata['anno']=allP_dict.keys()
allpdata['end']="*"


allpdatanode_seg = pd.DataFrame({'path': list(map(padict_letter.get, padict_letter.keys()))})
allpdatanode_seg['h']="P"
allpdatanode_seg['anno']=padict_letter.keys()
allpdatanode_seg['end']="*"





###unique haplotype
path_alsimlis=list(set(rowdictLETTERlist))
######################################generate micro gfafile##############################
##1.in:padict_letter NODEname+fa
node_name_dict={}
for nodename,path in padict_letter.items():
    symbols = re.findall(r'[+-]',path)
    numbers = re.findall(r'\d+', path)
    end = pd.DataFrame({'fa': list(map(fasta_dict.get, numbers))}, index=numbers)
    end['ori']=symbols[0:len(end)]
    fa_series = end.loc[end['ori'] == '-', 'fa']
    fa_series = fa_series.apply(lambda x: str(Seq(x).reverse_complement()))
    end.loc[end['ori'] == '-', 'fa'] = fa_series
    samseq=''.join(end['fa'])
    node_name_dict[nodename]=samseq


fasta_dict.update(node_name_dict)





##2.in:path_alsimlis






indicrev={}
for key,value in allnodes_dict.items():
    indicrev[key]=('|'.join(str.split(value,"|")[::-1]))

# inte=[]
# max_length_key = max(rowdictLETTER.keys(), key=lambda k: len(rowdictLETTER[k]))
# ben=rowdictLETTER[max_length_key]

# s=0
# e=s+2
# havall=[]
# while e<=len(ben):
#     letter='|'.join(ben[s:e])
#     vali=len(set([key for key, value in indic.items() if value.find(letter)!= -1]+[key for key, value in indicrev.items() if value.find(letter)!= -1]))
#     if vali==len(rowdict.keys()):
#         havall.append(letter)
#         e=e+1
#     else:
#         if len(havall)!=0:
#             inte.append(havall[-1])
#         s=e
#         e=s+2

# if len(havall)!=0:
#     inte.append(havall[-1])
# inte=set(inte)

# path_alsimlis=list(set(rowdictLETTERlist))
# intelis=path_alsimlis.copy()
# intenum=0
# for intenode in list(inte):
#     intelis = [s.replace(intenode, "INTE"+str(intenum)) for s in intelis]
#     intelis = [s.replace(("|".join(str.split(intenode,"|")[::-1])), "INTE"+str(intenum)) for s in intelis]
#     intenum=intenum+1

###继续整合
# intelisdict={}
# for i in range(len(intelis)):
#     intelisdict[i]=intelis[i]

# intelisdictrev={}
# for key,value in intelisdict.items():
#     intelisdictrev[key]=('|'.join(str.split(value,"|")[::-1]))

intelisdict=allnodes_dict
intelisdictrev=indicrev

path_alsimlis=list(set(uniqpath))
intelis=path_alsimlis.copy()
intelisdict={}
for i in range(len(intelis)):
    intelisdict[i]="|"+intelis[i]+"|"

intelisdictrev={}
for key,value in intelisdict.items():
    intelisdictrev[key]=('|'.join(str.split(value,"|")[::-1]))


expnode=[]
for m in intelis:
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


rowdictLETTERup={}
for index,value in rowdictLETTER.items():
    val=rowdictLETTER[index]
    rowdictLETTERup[index]="|".join([x for x in val if x not in ["+","-",""]])

alpath = pd.DataFrame(list(rowdictLETTERup.items()), columns=["Key", "Value"])
alpath.to_csv("out_alpath1"+region+".data",index=False,header=False,sep="\t")

intenum=0
for intenode in list(inte):
    for index,value in rowdictLETTERup.items():
        nname="INTE"+str(intenum)
        value=value.replace(intenode, nname) 
        rowdictLETTERup[index]=value
    intenum=intenum+1

alpath = pd.DataFrame(list(rowdictLETTERup.items()), columns=["Key", "Value"])
alpath.to_csv("out_alpath2"+region+".data",index=False,header=False,sep="\t")

other_indict_fa={}
other_indict_falis={}
intenum=0
for intenode in list(inte):
    nname="INTE"+str(intenum)
    other_indict_fa[nname]="".join(list(map(fasta_dict.get, str.split(intenode,"|"))))
    other_indict_fa[nname+"_rev"]=str(Seq(other_indict_fa[nname]).reverse_complement())
    intelis = [s.replace(intenode, nname) for s in intelis]
    intelis = [s.replace(("|".join(str.split(intenode,"|")[::-1])), nname) for s in intelis]
    other_indict_falis[nname]=intenode
    intenum=intenum+1

with open("out_INTE_node"+region+".json", "w", encoding="utf-8") as f:
    json.dump(other_indict_falis, f, ensure_ascii=False, indent=4)

fasta_dict.update(other_indict_fa)
allP_dictout,allL=inte_p(allP_dict,other_indict_falis)
allnode=(set("".join(allP_dictout.values()).replace("+",",").replace("-",",").split(","))-{""})

####整合path和inte:other_indict_falis allP_dict
# pal=inte_initnode(other_indict_falis,set(rowdictLETTERlist),padict_letter)
# pal=pd.concat([pal,allpdatanode_seg])
# filtered_pal = pal[pal['anno'].isin(allnode)]
# Pdata=filtered_pal[['h','anno','path','end']]
# Pdata.to_csv("P.data",index=False,header=False,sep="\t")

####生成inte的color
gene_color(allpdatanode_seg,other_indict_falis,0,allnode)
gene_color(allpdatanode_seg,other_indict_falis,1,allnode)
Pdata=allpdatanode_seg[['h','anno','path','end']]
Pdata.to_csv("NODE.P."+region+"data",index=False,header=False,sep="\t")


####生成gfa测试一下
allnode=set(allnode)
allL=pd.DataFrame(allL).drop_duplicates()
allpdata = pd.DataFrame({'path': list(map(allP_dictout.get, allP_dictout.keys()))})
allpdata['h']="P"
allpdata['anno']=allP_dictout.keys()
allpdata['end']="*"


Sdata = pd.DataFrame({'fa': list(map(fasta_dict.get, list(allnode)))})
# Sdata['fa']="A"
Sdata['node']=list(allnode)
Sdata['len']=Sdata['fa'].str.len()
Sdata['anno'] = Sdata['len'].apply(lambda x: f"LN:i:{x}")
Sdata['h']="S"
Sdata=Sdata[['h','node','fa','anno']]
Sdata.to_csv("S."+region+".data",index=False,header=False,sep="\t")

if len(allL)!=0:
    Ldata=allL.copy()
    Ldata.columns=['A','ori1','B','ori2']
    Ldata['h']="L"
    # Ldata['ori1']="+"
    # Ldata['ori2']="+"
    Ldata['anno']="0M"
    Ldata=Ldata[['h','A','ori1','B','ori2','anno']]
    Ldata['A_ori1'] = Ldata['A'].astype(str) + Ldata['ori1'].astype(str)
    Ldata['B_ori2'] = Ldata['B'].astype(str) + Ldata['ori2'].astype(str)
    Ldata['sorted_key'] = Ldata.apply(lambda row: tuple(sorted([row['A_ori1'], row['B_ori2']])), axis=1)
    Ldata['is_duplicate'] = Ldata.duplicated(subset=['sorted_key'], keep='first')
    Ldata = Ldata[~Ldata['is_duplicate']]
    Ldata=Ldata[['h','A','ori1','B','ori2','anno']]
    Ldata.to_csv("L."+region+".data",index=False,header=False,sep="\t")

Pdata=allpdata[['h','anno','path','end']]
# Pdata['path'] = Pdata['path'] + '+'
Pdata.to_csv("P."+region+".data",index=False,header=False,sep="\t")

print(padict_letter)
out=dict_keys_jaccard_similarity(padict_letter)

pd.DataFrame(out).to_csv("similarity_node"+region+".bed",index=False,header=False,sep="\t")
