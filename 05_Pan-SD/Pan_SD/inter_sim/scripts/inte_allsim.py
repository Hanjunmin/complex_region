import pandas as pd
import argparse
import pandas as pd


df = pd.read_csv("reg.sdsim.2k",sep="\t", index_col=0)
x1 = df[(df['1'] == ".") & (df['2'] != ".")]
x2=df[(df['1'] != ".") & (df['2'] == ".")]
x3=df[(df['1'] == ".") & (df['2'] == ".")]
x4=df[(df['1'] != ".") & (df['2'] != ".")]
if len(x1)!=0:
	x1['sim']=x1['2']

if len(x2)!=0:
	x2['sim']=x2['1']

if len(x3)!=0:
	x3['sim']="."

if len(x4)!=0:
	x4['sim']=x4['1'].astype("float")+x4['2'].astype("float")

x4.loc[x4['sim']>1,'sim']=1 ##直接截断
end=pd.concat([x1,x2,x3,x4])


# dictal={}
# for win in set(list(end['5'])):
# 	print(win)
# 	subdf=end[end['5']==win]
# 	key_value_dict = dict(zip(subdf.iloc[:, 0], subdf.iloc[:, 6]))
# 	dictal[win]=key_value_dict
dictal = end.groupby('5').apply(
    lambda g: dict(zip(g.iloc[:, 0], g.iloc[:, 6]))
).to_dict()

end=pd.DataFrame(dictal)
simdf=end.transpose()
simdf['win']=simdf.index

# ano1=pd.read_csv('/home/jmhan/PANSDEND/000_SD_sim/win/'+CHR+'.simn_ext',sep="\t", index_col=0)
# ano2=pd.read_csv('/home/jmhan/PANSDEND/000_SD_sim/win/'+CHR+'.simn_long',sep="\t", index_col=0)
# Z=pd.concat([ano1,ano2])
# Z['win']=Z.index
# Z_filtered = Z[~Z['win'].isin(simdf['win'])]


ano3=pd.read_csv('reg.out',sep="\t", index_col=0)
ano3['win']=ano3.index
ano3_filtered =ano3[~ano3['win'].isin(simdf['win'])]




END = pd.concat([simdf, ano3_filtered]) 
END.to_csv('sim.out',sep='\t')