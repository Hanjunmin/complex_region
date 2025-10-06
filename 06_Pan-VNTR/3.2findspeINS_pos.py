import pysam
import pandas as pd
import pickle
from collections import defaultdict
import glob
import argparse
import pysam
import sys
sys.path.append('./project/ENDVNTR')
from vntr_fun import *
from Bio.Seq import Seq
import warnings
import subprocess

warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser(description='OTHER-SDcoordi')
parser.add_argument('--sam', help='sample')
args = parser.parse_args()


def generate_motif_pos(iterid):
    motiflen=len(cnmotif)
    st=stal-iterid*motiflen-20
    en=enal+iterid*motiflen+20
    if st==0:
        st=st+1
    if st<0:
        st=1
    if en>maxlength:
        en=maxlength
    T2TSeq = fa.fetch(region=subchr+":"+str(st)+"-"+str(en))
    ref=T2TSeq.upper()
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
    return storelistn,st

def refinecor(num):
    if num<1:
        num=1
    if num>maxlength:
        num=maxlength
    return num
    



# seq="GGGTTGGTCCTGCCCCTCAGAGGCAGAAAGGGCAGAAGGGGCAGAAGAG*******GTGGAGAAGGTGGAAGCGGGGCAGAAGAG*******GTGGAGAAGGTGGAAGGGAGGCAGGAGAGGCAGGCACATTCAGTGTAAC"
# seq=seq.replace("*", "")
# motif="GGGGCAGAAGAGGTGGAGAAGGTGGAAG"

# sam="C001-CHA-E01-Mat"        
# subchr="C001-CHA-E01_Mat_chr10"
# s=124354061     
# e=124354077



sam=args.sam   #"NA19650_hap2"
print(sam)
#sam="GRCh38"
FA="./allfa/"+sam+".fa"
fa = pysam.FastaFile(FA)
REFfa = pysam.FastaFile("./allfa/CHM13v2.fa")

with open(sam+".out", 'r') as f:
    lines = f.readlines()

from collections import Counter

id=0
x=[]
alldflis=[]
for line in lines[1:]:
    line_list = line.strip().split('\t')
    subchr=line_list[14]
    subchrs=int(line_list[15])
    maxlength=fa.get_reference_length(subchr)
    posa=[]
    for subl in line_list[13].split("@"):
        posa.append((int(subl.split("_")[0]),int(subl.split("_")[1])))
    have="FALSE"
    for possam in posa:
        s=possam[0]+subchrs
        e=possam[1]+subchrs
        insseq=line_list[7].split("*******")[1]
        seq=line_list[7].replace("N****", "")
        seq=seq.replace("*****N", "")
        seq=seq.replace("*", "")
        motif=line_list[4]
        first=refinecor(s-len(seq))
        second=refinecor(e+len(seq))
        T2TSeq = fa.fetch(region=subchr+":"+str(first)+"-"+str(second))
        T2TSeq=T2TSeq.upper()
        if T2TSeq.find(seq) != -1:
            seq=seq
            have="TRUE"
            break
        # elif T2TSeq.find(str(Seq(seq).reverse_complement()))!=-1:
        #     seq=str(Seq(seq).reverse_complement())
    if have=="FALSE":
        id=id+1
        x.append(line)
        # print(sam)
        # print(line)
        # print("NO")
        break
    stal=first+T2TSeq.find(seq)
    enal=stal+len(seq)
    T2TSeq = fa.fetch(region=subchr+":"+str(stal)+"-"+str(enal))
    cnmax=0
    cnmotif="a"
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
            cnmax=cn1
            cnmotif=motif
    if alllenlist2!=[]:
        if cn2>cnmax:
            cnmax=cn2
            cnmotif=str(Seq(motif).reverse_complement())
    if cnmax==0:
        poslis=[line_list[0],line_list[1],line_list[2],0,motif]
        alldflis.append(poslis)
        continue
    else:
        initcn=cnmax  
        nextcn=initcn+1
        initid=1
    while nextcn>initcn:
        [storelist,st]=generate_motif_pos(initid)
        out1=process_vamos(storelist)
        nextcn=out1[0][0]
        initid=initid+2
        if nextcn>initcn:
            initcn=nextcn
    molis=[out1[0][4]] if not isinstance(out1[0][4], list) else set(out1[0][4])
    positions = []
    for start, end in out1[0][3]:
        positions.extend(range(start, end + 1))  # 展开区间为单个位置
    position_counts = Counter(positions)
    count_one_positions = [pos for pos, count in position_counts.items() if count == 1]
    simindex=len(count_one_positions)/((out1[0][2]-out1[0][1])+1)
    poslis=[line_list[0],line_list[1],line_list[2],subchr,nextcn,sam,len(list(molis)[0]),st+out1[0][1],st+out1[0][2],",".join(molis),stal,enal,seq.find(insseq),seq.find(insseq)+len(insseq),line_list[6],line_list[9],simindex] ###最后四列为样本整个seq的区域，和ins在seq中的位置——seq为找到的ins和左右anchor的区域,第7列第八列为motif回去比对左右扩后得到的坐标 最后生成了一个sim
    alldflis.append(poslis)

        
dfal=pd.DataFrame(alldflis)
sammaxlength=[fa.get_reference_length(subchr) for subchr in dfal[3]]
refmaxlength=[REFfa.get_reference_length(subchr) for subchr in dfal[0]]
dfal['sammaxlen']=sammaxlength
dfal['REFmaxlen']=refmaxlength
dfal['len']=dfal[11]-dfal[10]
dfal[12]=dfal[10]+dfal[12]
dfal[13]=dfal[10]+dfal[13]
dfal.columns=['REFCHR','REFPOS','TYPE','SAM_CHR','CN','SAM','MOTIF_LEN','RP_samS','RP_samE','MOTIF','Seq_samS','Seq_samE','INS_samS','INS_samE','sampath','cluster','motifsim','sammaxlen','REFmaxlen','Seqlen']
'''
RP_samS:真正VNTR的位置
RP_samE:真正VNTR的位置
Seq_samS:Seq(graph上的ins+左右REF anchor)的位置
Seq_samE:Seq(graph上的ins+左右REF anchor)的位置
INS_samS:Seq中ins的位置
INS_samE:Seq中ins的位置

'''
dfal['len']=dfal['INS_samE']-dfal['INS_samS']
####分为两种情况 一个是左右扩去找motif对应的区域的时候扩到seq之外了(expandye)  一个是没有扩到外面（expandno）
expandye=dfal[(dfal['RP_samS']<dfal['Seq_samS']) | (dfal['RP_samE']>dfal['Seq_samE'])]
expandye=expandye[(expandye['len']==5000)] ###找到500且左右repeats扩出去了的
if len(expandye.index)!=0:
    expandno = dfal.drop(expandye.index)
else:
    expandno=dfal

###先看一下expandno
###1.比对验证是否是ins
def minimap_line(linen,foldn):
    sams=line['Seq_samS']-foldn*line['Seqlen']
    same=line['Seq_samE']+foldn*line['Seqlen']
    refs=int(line['REFPOS'])-foldn*line['Seqlen']
    refe=int(line['REFPOS'])+foldn*line['Seqlen']
    samsn=max(1,sams)
    samen=min(same,line['sammaxlen'])
    refsn=max(1,refs)
    refen=min(refe,line['REFmaxlen'])
    sam=line['SAM']
    samchr=line['SAM_CHR']
    refchr=line['REFCHR']
    A_cmd = f"samtools faidx ./allfa/CHM13v2.fa {refchr}:{refsn}-{refen} > {sam}ref.fa"
    # print(A_cmd)
    sam_cmd = f"samtools faidx ./allfa/{sam}.fa {samchr}:{samsn}-{samen} > {sam}.fa"
    # print(sam_cmd)
    subprocess.run(A_cmd, shell=True, check=True)
    subprocess.run(sam_cmd, shell=True, check=True)
    minimap_cmd = f"minimap2 -cx asm5 {sam}ref.fa {sam}.fa >{sam}.paf"
    result = subprocess.run(minimap_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = result.stdout.decode('utf-8')
    with open(sam+".paf", 'r') as file:
        lines = file.readlines()
    if len(lines)==0:
        return samsn,samen,refsn,refen,"NOALIGN"
    else:
        return samsn,samen,refsn,refen,"ALIGN"
    


def split_into_chunks(nums):
    s = pd.Series(nums)
    diff = s.diff().ne(1).cumsum()
    return [(group.min(), group.max()) for _, group in s.groupby(diff)]


def get_paf(start,linen):
    with open(sam+'.paf', 'r') as file:
        lines = file.readlines()
    if len(lines)==0:
        return 0,"NO"
    paflist=[]
    haveset=[]
    for pafline in lines:
        line_list = pafline.strip().split('\t')
        cg_items = [item for item in line_list if item.startswith('cg:')]
        subline=line_list[0:9]
        subline.extend(cg_items)
        matches = re.findall(r'\d+I', cg_items[0])
        if matches:
            numbers = [int(match[:-1]) for match in matches]  # Remove 'I' and convert to integer
            max_number = max(numbers)
        else:
            max_number =0
        subline.extend([max_number])
        paflist.append(subline)
        haveset.extend(range(int(line_list[2]),int(line_list[3])))
    pafdf=pd.DataFrame(paflist)
    a=min(pafdf[2].astype(int))
    b=max(pafdf[3].astype(int))
    subnum=set(range(a,b))-set(haveset) # out=split_into_chunks(list(subnum))
    # print(pafdf)
    inslenn=max(max_number,len(subnum))
    if linen['INS_samS']>=a+start & linen['INS_samE']<=b+start: ###这里很粗略，因为可能周围有别的INS，但是我们也不知道
        pafinsin="TRUE"
    else:
        pafinsin="FALSE"
    ratio=min(inslenn,linen['insbef'])/max(inslenn,linen['insbef'])
    return ratio,pafinsin

expandno['INSinPAF']="NO"
expandno['ALIGN']="NO"
expandno['INSlen']=0
expandno['insbef']=expandno['INS_samE']-expandno['INS_samS']
for index,line in expandno.iterrows():
    samsn,samen,refsn,refen,ALIGNTYPE=minimap_line(line,1)
    # print(line['REFCHR']+":"+str(refsn)+"-"+str(refen))
    ratio,alipos=get_paf(samsn,line)
    expandno.loc[index, 'INSinPAF'] = alipos
    expandno.loc[index, 'ALIGN'] = ALIGNTYPE
    expandno.loc[index, 'INSratio'] = ratio
    if ratio<0.1:
        fold=1
        while ratio<0.1:
            if fold>6:
                break
            fold=fold+1
            samsn,samen,refsn,refen,ALIGNTYPE=minimap_line(line,fold)
            ratio,alipos=get_paf(samsn,line)
        expandno.loc[index, 'INSratio'] = ratio
        expandno.loc[index, 'INSinPAF'] = alipos
        expandno.loc[index, 'ALIGN'] =ALIGNTYPE

###expandno 保留ratio>0.2的
if len(expandno)!=0:
    savesmvntr=expandno[expandno['INSratio']>=0.2]
    savesmvntr.to_csv(sam+'smvntr.out', index=False,sep="\t")

####下面为expandyes
expandye['INSinPAF']="NO"
expandye['ALIGN']="NO"
expandye['INSlen']=0
expandye['insbef']=expandye['INS_samE']-expandye['INS_samS']
for index,line in expandye.iterrows():
    samsn,samen,refsn,refen,ALIGNTYPE=minimap_line(line,1)
    # print(line['REFCHR']+":"+str(refsn)+"-"+str(refen))
    ratio,alipos=get_paf(samsn,line)
    expandye.loc[index, 'INSinPAF'] = alipos
    expandye.loc[index, 'ALIGN'] = ALIGNTYPE
    expandye.loc[index, 'INSratio'] = ratio
    if ratio<0.1:
        fold=1
        while ratio<0.1:
            if fold>3:
                break
            fold=fold+1
            samsn,samen,refsn,refen,ALIGNTYPE=minimap_line(line,fold)
            ratio,alipos=get_paf(samsn,line)
        expandye.loc[index, 'INSratio'] = ratio
        expandye.loc[index, 'INSinPAF'] = alipos
        expandye.loc[index, 'ALIGN'] =ALIGNTYPE

if len(expandye)!=0:
    savelarvntr=expandye[expandye['INSratio']>=0.1]
    savelarvntr.to_csv(sam+'larvntr.refine', index=False,sep="\t")



