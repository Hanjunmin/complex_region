import pandas as pd
import argparse
import os 
parser = argparse.ArgumentParser(description='constructdict')
parser.add_argument('--pair', help='PAIR PATH')
parser.add_argument('--sam', help='sample name')
parser.add_argument('--scp', help='R script')

args = parser.parse_args()


# file_path = '/home/jmhan/project/APG/github_final/PanSD_endtest/run/03_nonref/02_ITER/2_SD_GENE/tranall.uniq'
file_path=args.pair
pepair=pd.read_csv(file_path,sep="\t",header=None)
pepair.columns=['A', 'B']


sam=args.sam
allpair1=pd.read_csv(sam+".pos_sim",sep="\t",header=0)
mask = allpair1['A.pos'] == 'NOSAM'
allpair1.loc[mask, 'A.pos'] = 'not_align'
allpair1.loc[mask, 'A.sim'] = 0

mask = allpair1['B.pos'] == 'NOSAM'
allpair1.loc[mask, 'B.pos'] = 'not_align'
allpair1.loc[mask, 'B.sim'] = 0



allp1=pd.merge(pepair,allpair1, how='left', on=["A","B"]).drop_duplicates()
allp1.columns=["chr1_init","chr2_init","chr1_source","chr1_init.sim","chr2_source","chr2_init.sim"]
duplicate_mask = allp1.duplicated(subset=['chr1_init', 'chr2_init'], keep=False)

# 然后修改这些重复行的sim列
allp1.loc[duplicate_mask, 'chr2_init.sim'] = '.'
allp1.loc[duplicate_mask, 'chr1_init.sim'] = '.'
allp1.loc[duplicate_mask, 'chr1_source'] = '.'
allp1.loc[duplicate_mask, 'chr2_source'] = '.'

allp1=allp1.drop_duplicates()

###去重
file_size = os.path.getsize(sam+".eebed.rev.out")

if file_size == 0:
	END=allp1
else:
	ENDpari=pd.read_csv(sam+".eebed.rev.out",sep="\t",header=None)
	ENDparisav=ENDpari.iloc[:,0:12].drop_duplicates()
	ENDparisav.columns=["chr1","start1","end1","chr2","start2","end2","match_base","mismatch_base","orient","identity","align_len","matchindel"]
	ENDparinit=pd.read_csv(sam+".all1.end.statistics", skiprows=1,sep="\t",header=None)
	ENDparinit.columns=["chr1","start1","end1","chr2","start2","end2","match_base","mismatch_base","orient","identity","align_len","matchindel","chr2_source","chr1_source","sim1","sim2"]
	ENDparinit=ENDparinit.iloc[:,0:14].drop_duplicates()
	allp2=pd.merge(ENDparisav,ENDparinit, how='left', on=["chr1","start1","end1","chr2","start2","end2","match_base","mismatch_base","orient","identity","align_len","matchindel"])
	allp2.to_csv(sam+".eebed.rev.out.ad", index=False,sep="\t",header=False)
	colname=allp2.columns
	import subprocess
	result = subprocess.run(
	        ["Rscript", args.scp, sam+".eebed.rev.out.ad", sam+".eebed.rev.out.out"],
	        capture_output=True,
	        text=True,
	        check=True
	    )
	allp2=pd.read_csv(sam+".eebed.rev.out.out",sep="\t",header=None)
	allp2.columns=colname
	END=pd.merge(allp1,allp2,how='left', on=["chr2_source","chr1_source"])


END = END.fillna('.')
END['SAM']=sam
END.to_csv(sam+".ITER.SD.al", index=False,sep="\t")

