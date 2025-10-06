import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import numpy as np
import skfuzzy as fuzz
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description='RUNMODEL')
parser.add_argument('--m', help='modelname(KM OR FCC)')
parser.add_argument('--c', help='cluster number')
parser.add_argument('--b', help='bubbleType(1 or 2)')
args = parser.parse_args()
model=args.m
kc=int(args.c)
bubty=int(args.b)

model="FCC"
kc=5
bubty=2

initdata2 = pd.read_csv("./COMPLEXFEATURE/05_split/allfeaturen",header=0)  ##allfeaturen是后来改了vscore的结果
initdata2.columns=['chr','start','end','SDsim','vntr','Vscore','entropy','Vscoren','simple','super','insertion','anno_vscore_na','anno_simple_na']

refinebubble = pd.read_csv("./COMPLEXFEATURE/01_bubble/para_gene_path2/end.txt",header=None,sep="\t")
if bubty==1:
    refinebubble.columns=['chr','start','end','simple','super','insertion','poly','sim','simple_change','super_change','insertion_change','simple_change1','super_change1','insertion_change1','z','x','y']
else:
    refinebubble.columns=['chr','start','end','simple','super','insertion','poly','sim','simple_change1','super_change1','insertion_change1','simple_change','super_change','insertion_change','z','x','y']

refinebubble['start']=refinebubble['start']-1
initdatac=pd.merge(initdata2, refinebubble, how='left', on=['chr','start','end','simple','super','insertion'])
initdatac.loc[initdatac['simple_change'].notna(), 'simple'] = initdatac['simple_change']
initdatac.loc[initdatac['super_change'].notna(), 'super'] = initdatac['super_change']
initdatac.loc[initdatac['insertion_change'].notna(), 'insertion'] = initdatac['insertion_change']

initdatac=initdatac.loc[:,['chr','start','end','SDsim','vntr','Vscore','entropy','Vscoren','simple','super','insertion']]
# APG=initdata[initdata['chr']=="chr1"]
APG2=initdatac.copy()
APG2['Vscore']=initdatac['Vscoren']*initdatac['entropy']
APG2['simplen']=0.6*APG2['simple']+0.1*APG2['insertion']+0.3*APG2['super']
APG2['alln'] = np.log1p(APG2['simplen']) 

APG2['simplen'] = np.sqrt(APG2['simple']) 
APG2['supern'] =np.sqrt(APG2['super'])  
APG2['insertionn'] =np.sqrt(APG2['insertion'])  

# APG2=APG2[APG2['simple']<20000]
# APG2.loc[APG2['vntr'] < 0.05, 'vntr'] = 0
# APG2.loc[APG2['simple'] <30 , 'insertion'] = 0
# APG2.loc[APG2['simple'] <30 , 'super'] = 0
# APG2.loc[APG2['simple'] <30 , 'simple'] = 0
# APG2.loc[APG2['Vscore'] <0.14 , 'Vscore'] = 0
APGcopy2=APG2.copy()
APGbef2=APGcopy2.iloc[:,[0,1,2]]
APGin2=APGcopy2[['SDsim','simplen','supern','insertionn','Vscore','vntr']]
# APGin2=APGcopy2[['SDsim','simplen','Vscore','vntr']]
# APGin2.loc[APGin2['SDsim'] < 0.2, 'SDsim'] = 0



#####kmeans model
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import skfuzzy as fuzz
from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics import accuracy_score
import sklearn.metrics
import sklearn.metrics as metrics
import pandas as pd
if model=="KM":
    X = APGin2
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)  # 标准化处理
    kmeans = KMeans(n_clusters=kc, n_init=20)
    y_pred = kmeans.fit_predict(X_scaled)
    labels =y_pred
    APGcopy2['cluster']=labels
    # complex_regchr.index=pd.DataFrame(u.T).index
    # allnew=pd.concat([complex_regchr[['chr','start','end','SDsim','simple','super','insertion','Vscore','vntr','cluster']],pd.DataFrame(u.T),],axis=1)
    plt.figure(figsize=(14, 8))
    for i, feature in enumerate(['Vscore', 'SDsim', 'simple', 'super', 'insertion','vntr']):
        plt.subplot(2, 3, i+1)
        sns.boxplot(x='cluster', y=feature, data=APGcopy2)
        plt.title(f'Distribution of {feature} by cluster')
    plt.tight_layout()
    plt.savefig('./COMPLEXFEATURE/08_model/'+model+"_clus"+str(kc)+"@"+str(bubty)+".png", format='png', bbox_inches='tight')
    plt.show()
    end=APGcopy2.loc[:,['chr','start','end','SDsim','vntr','Vscore','simple','super','insertion','cluster']]
    end.to_csv('./COMPLEXFEATURE/08_model/'+model+"_clus"+str(kc)+"@"+str(bubty)+".csv")

SSE = [] 
for k in range(1, 9):
    estimator = KMeans(n_clusters=k) 
    estimator.fit(X_scaled)
    SSE.append(estimator.inertia_)




#####fcc model
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import skfuzzy as fuzz
from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics import accuracy_score
import sklearn.metrics
import sklearn.metrics as metrics
import pandas as pd
if model=="FCC":
    X = APGin2
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)  # 标准化处理
    # kmeans = KMeans(n_clusters=5, n_init=20)
    # y_pred = kmeans.fit_predict(X_scaled)
    # labels =y_pred
    cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(
        X_scaled.T,kc, 1.5, error=0.05, maxiter=1000, init=None)
    cluster_membership = np.argmax(u, axis=0)
    labels = cluster_membership
    APGcopy2['cluster']=labels
    end=pd.concat([APGcopy2[['chr','start','end','SDsim','vntr','Vscore','simple','super','insertion','cluster']],pd.DataFrame(u.T),],axis=1)
    # complex_regchr.index=pd.DataFrame(u.T).index
    # allnew=pd.concat([complex_regchr[['chr','start','end','SDsim','simple','super','insertion','Vscore','vntr','cluster']],pd.DataFrame(u.T),],axis=1)
    plt.figure(figsize=(14, 8))
    for i, feature in enumerate(['Vscore', 'SDsim', 'simple', 'super', 'insertion','vntr']):
        plt.subplot(2, 3, i+1)
        sns.boxplot(x='cluster', y=feature, data=APGcopy2)
        plt.title(f'Distribution of {feature} by cluster')
    plt.tight_layout()
    plt.savefig('./COMPLEXFEATURE/08_model/'+model+"_clus"+str(kc)+"@"+str(bubty)+".png", format='png', bbox_inches='tight')
    plt.show()
    end.to_csv('./COMPLEXFEATURE/08_model/'+model+"_clus"+str(kc)+"@"+str(bubty)+".csv")



import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture

custom_palette = ['#A2408B','#806C9D',  '#CAAFC2',  '#6D96BF','#2A5192', '#BBB9B9']
  

plt.figure(figsize=(10, 4))  
for i, feature in enumerate(['SDsim', 'vntr','simple', 'super', 'insertion','Vscore' ]):
    plt.subplot(1, 7, i+1) 
    sns.violinplot(x=feature, y='cluster', data=end, inner='box', palette=custom_palette,orient='h',width=0.8)
    plt.xlabel(feature)  
    if i > 0:
        plt.ylabel('')
        plt.yticks([])
    if i == 0:
        plt.ylabel('Cluster') 
        plt.yticks(range(len(set(end['cluster']))), sorted(set(end['cluster']))) 
plt.tight_layout()
plt.savefig('./COMPLEXFEATURE/08_model/ig1_deepclustermodel.violin.pdf', format='pdf', bbox_inches='tight')
plt.show()
