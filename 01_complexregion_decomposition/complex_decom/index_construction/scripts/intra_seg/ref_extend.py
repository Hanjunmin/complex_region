import pandas as pd
import pickle
from collections import defaultdict
import glob
import argparse
import pysam
import sys
from collections import Counter
from Bio.Seq import Seq
import warnings
warnings.filterwarnings('ignore')







def generate_motif_pos(iterid,trans):
    motiflen=len(motifdic[row[0]][0])
    st=stal-iterid*motiflen-20
    en=enal+iterid*motiflen+20
    if st==0:
        st=st+1
    if st<0:
        st=1
    if en>maxlength:
        en=maxlength
    T2TSeq = fa.fetch(region=row['chr']+":"+str(st)+"-"+str(en))
    ref=T2TSeq.upper()
    if trans==2: ##chan2
        ref=str(Seq(T2TSeq).reverse_complement())
    storelistn=[]
    motif=cnmotif
    resdict,interval, motiflist, alllenlist = defaultdict(list),(0, len(ref)), [motif], []
    expand_motif(resdict,motif, ref, interval,0.2,alllenlist)
    if alllenlist!=[]:
        cn=cal_cn1(alllenlist,len(motif),0.3)
        cn=list(cn)
        motiflist=[]
        for mo in cn[3]:
            motiflist.append(ref[mo[0]:(mo[1]+1)])
        cn.append(motiflist)
        storelistn.append(cn)   
    return storelistn,st,en-st ##chan1

parser = argparse.ArgumentParser("integrate pkl file")
# parser.add_argument("-init", "--init", type=str, required=True, help="the init cn file")
parser.add_argument("-i", "--i", type=str, required=True, help="the input file")
parser.add_argument("-o", "--o", type=str, required=True, help="the output file")
parser.add_argument('--script', help='script of VNTR calculation')
parser.add_argument('--fapath', help='VNTR')
args = parser.parse_args()



sys.path.append(args.script) #'/home/jmhan/project/APG/github_final/01_1p36/index_construction/script/intra_seg'
from vntr_fun import *


# initcn=args.init
file_name=args.i
# file_name="0pos.pkl"
initcn="can.VNTR.pos"
T2Tpos = pd.read_csv(initcn,sep="\t",header=None)
T2Tpos.columns=['CHR','S','E','CN','MOTIF']
T2Tpos['inde'] = T2Tpos['CHR'] + "_" + T2Tpos['S'].astype(str) + "_" + T2Tpos['E'].astype(str)
motifdic=T2Tpos.set_index('inde')['MOTIF'].str.split(",").to_dict()




# def process_pkl(file_name):

print(file_name)
with open(file_name, 'rb') as f:
    results = pickle.load(f)  # 从文件中加载结果

alldflis=[]
if len(results)!=0:
    xx=[]
    for key,value in results.items():
        for samal,pos in value.items():
            xx.append([key,samal.split("@")[0],samal.split("@")[2],samal.split("@")[3],samal.split("@")[4],pos[0],pos[1]])
    xxdf=pd.DataFrame(xx)
    xxdf['inits'] = pd.to_numeric(xxdf[0].str.split("_").str[1])
    xxdf['inite'] = pd.to_numeric(xxdf[0].str.split("_").str[2])
    xxdf['ns'] = xxdf[5]
    xxdf['ne'] = xxdf[6]
    xxdf['ratio'] = abs(xxdf['inite'] - xxdf['inits']) / np.maximum(
    abs(xxdf['inite'] - xxdf['inits']), 
    abs(xxdf['ne'] - xxdf['ns']))
    xxdf['chr']=xxdf[2]+"@"+xxdf[3]
    # delregion=xxdf[xxdf['ratio']<0.001]  ##modify 这里注释了很多 感觉不需要del 因为是VNTR 直接后面加一个如果序列整体长度大于500000就continue吧
    # for index,row in delregion.iterrows():
    #     molis=[motifdic[row[0]]] if not isinstance(motifdic[row[0]], list) else motifdic[row[0]]
    #     poslis=[row[0],"LARRANGE",row[1],row[2],0,0,",".join(molis)]
    #     alldflis.append(poslis)
    # xxdf=xxdf[xxdf['ratio']>=0.01]
    for sam in set(xxdf[1]):
        print(sam)
        FA=args.fapath+sam+".fa"
        fa = pysam.FastaFile(FA)
        xxdfsam=xxdf[xxdf[1]==sam] 
        for index,row in xxdfsam.iterrows():
            # print(index)
            maxlength=fa.get_reference_length(row['chr'])
            motiflen=len(motifdic[row[0]][0])
            stal=min(row['ns'],row['ne'])
            if stal==0:
                stal=1
            enal=max(row['ns'],row['ne'])
            ####确定初始motif和对应的拷贝数
            T2TSeq = fa.fetch(region=row['chr']+":"+str(stal)+"-"+str(enal))
            T2TSeq=T2TSeq.upper()
            cnmax=0
            cnmotif="a"
            for motif in motifdic[row[0]]:
                ref=T2TSeq
                resdict,interval, motiflist, alllenlist1 = defaultdict(list),(0, len(ref)), [motif], []
                expand_motif(resdict,motif, ref, interval,0.2,alllenlist1)
                samseqrev=str(Seq(T2TSeq).reverse_complement())
                ref=samseqrev
                resdict,interval, motiflist, alllenlist2 = defaultdict(list),(0, len(ref)), [motif], []
                expand_motif(resdict,motif, ref, interval,0.2,alllenlist2)
                if alllenlist1!=[]:
                    enall=cal_cn1(alllenlist1,len(motif),0.3)
                    cn1=enall[0]
                if alllenlist2!=[]:
                    enall=cal_cn1(alllenlist2,len(motif),0.3)
                    cn2=enall[0]
                if alllenlist1!=[]:
                    if cn1>cnmax:
                        trans=1
                        cnmax=cn1
                        cnmotif=motif
                        # vntridcopy_number[sam]=cn1
                if alllenlist2!=[]:
                    if cn2>cnmax:
                        trans=2
                        cnmax=cn2
                        #cnmotif=str(Seq(motif).reverse_complement())
                        cnmotif=motif
            if cnmax==0:
                poslis=[row[0],0,sam,row['chr'],row['ne'],row['ns'],motif,-1]
            else:
                ini=0          
                for cnmotif in motifdic[row[0]]: ###相当于找到最concensus的motif序列然后进行注释
                    initcn=cnmax  
                    outmax="a"
                    nextcn=initcn+1
                    initid=1
                    lastcn=0
                    conti=0
                    while nextcn>initcn & conti<10:
                        if initid!=1:
                            outmax=out1
                        [storelist,st,lena]=generate_motif_pos(initid,trans) ##chan3
                        if len(storelist)==0:
                            break
                        out1=process_vamos(storelist)
                        out1.append(lena)
                        # print(out1)
                        nextcn=out1[0][0]
                        initid=initid+2
                        if nextcn==lastcn  or nextcn<lastcn:
                            conti=conti+1
                        if nextcn>lastcn:
                            cncopy=out1
                        if conti>10:
                            outmax=cncopy
                            break
                        lastcn=nextcn
                    if initid==3:
                        outmax=out1
                    # print(lastcn)
                    if lastcn>ini:
                        ini=lastcn
                        cnn=outmax
                outmax=cnn
                # print(cnn)
                molis=[outmax[0][4]] if not isinstance(outmax[0][4], list) else set(outmax[0][4])
                if outmax[0][2]-outmax[0][1]!=0:
                    positions = []
                    for start, end in outmax[0][3]:
                        positions.extend(range(start, end + 1))  # 展开区间为单个位置
                    position_counts = Counter(positions)
                    count_one_positions = [pos for pos, count in position_counts.items() if count == 1]
                    simindex=len(count_one_positions)/((outmax[0][2]-outmax[0][1])+1)
                    #simindex=(sum(end - start for (start, end) in outmax[0][3]))/(outmax[0][2]-outmax[0][1])
                else:
                    simindex=-1
                if trans==2:
                    poslis=[row[0],outmax[0][0],sam,row['chr'],st+outmax[-1]-outmax[0][2],st+outmax[-1]-outmax[0][1],",".join(molis),simindex]
                else:
                    poslis=[row[0],outmax[0][0],sam,row['chr'],st+outmax[0][1],st+outmax[0][2],",".join(molis),simindex]
            alldflis.append(poslis)
    alldflisdf=pd.DataFrame(alldflis)
    alldflisdf.to_csv(args.o, index=True, sep="\t")
else:
    empty_df = pd.DataFrame()
    empty_df.to_csv(args.o, index=False, sep="\t")
