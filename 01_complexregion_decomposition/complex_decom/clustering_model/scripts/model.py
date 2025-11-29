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
from time import time
import numpy as np
import keras.backend as K
from tensorflow.keras.layers import Layer, InputSpec #from keras.engine.topology import Layer, InputSpec
from keras.layers import Dense, Input
from keras.models import Model
from keras.optimizers import SGD
from keras import callbacks
from keras.initializers import VarianceScaling
from sklearn.cluster import KMeans
import tensorflow as tf
import argparse
parser = argparse.ArgumentParser(description='RUNMODEL')

parser.add_argument('--dir', help='base dir')
parser.add_argument('--model', help='model parameters path')
args = parser.parse_args()

outputdir=args.dir
Modeldir=args.model


win=pd.read_csv(outputdir+"/01_preproc/cor_window.bed",sep="\t",header=None)  ##allfeaturen是后来改了vscore的结果
win.columns=['chr','start','end']


###Seq_div
entropy=pd.read_csv(outputdir+"02_index_construction/seq_div/entropy_final.txt",sep="\t",header=None)
entropy.columns=['chr','start','end',"entropy"]
FEATUREIN=entropy
simial=pd.read_csv(outputdir+"02_index_construction/seq_div/simal",sep="\t",header=0)
simial.columns=['chr','start','end',"win_sim"]
FEATUREIN=pd.merge(FEATUREIN, simial, how='left', on=['chr','start','end'])

###TE
TEvar=pd.read_csv(outputdir+"02_index_construction/TE_div/TEEND",sep="\t",header=0)
TEvar=TEvar[['chr','start','end','trimgini']]
TEvar.columns=['chr','start','end',"TEgini"]
TEvar['TEgini'] = TEvar['TEgini'].fillna(0)
FEATUREIN=pd.merge(FEATUREIN, TEvar, how='left', on=['chr','start','end'])

###inter_seg
interseg=pd.read_csv(outputdir+"02_index_construction/inter_seg/reg.simout",sep="\t",header=0)
interseg[['chr', 'start', 'end']] = interseg['win'].str.split(r'[:-]', expand=True)
interseg=interseg[['chr','start','end','sdrow_max','sdtrimsd']]
interseg.columns=['chr','start','end',"sd_max",'sd_var']
interseg['start']=interseg['start'].astype("int")
interseg['end']=interseg['end'].astype("int")


FEATUREIN=pd.merge(FEATUREIN, interseg, how='left', on=['chr','start','end'])


###intra_seg
intraseg= pd.read_csv(outputdir+"02_index_construction/intra_seg/intra_seg.fea",header=0,sep="\t")  ##allfeaturen是后来改了vscore的结果
intraseg.loc[intraseg['vntr_cov']>0,"vntr_cov"]=1
intraseg=intraseg[['chr','start','end','vntr_cov','vntr_gini_max']]
FEATUREIN=pd.merge(FEATUREIN, intraseg, how='left', on=['chr','start','end'])


###bubble
FEATUREIN['sim_bubble']=0.01
FEATUREIN['com_bubble']=0.01
refinebubble = pd.read_csv(outputdir+"02_index_construction/bub/bub.txt",header=None,sep="\t")
refinebubble.columns=['chr','start','end','chrregion','startregion','endregion','poly','polyratio','simple','pathratio','lenration','pos','polyrefine','H','V4new','V5new','norm_bub','bub']
refinebubble['sim_bub']=refinebubble['V5new']+refinebubble['H']
refinebubble['com_bub']=refinebubble['V4new']+refinebubble['H']
FEATUREIN=pd.merge(FEATUREIN, refinebubble[['chr','start','end','sim_bub','com_bub']], how='left', on=['chr','start','end'])
FEATUREIN.loc[FEATUREIN['sim_bub'].notna(), 'sim_bubble'] = FEATUREIN['sim_bub']
FEATUREIN.loc[FEATUREIN['com_bub'].notna(), 'com_bubble'] = FEATUREIN['com_bub']





FEATUREIN=FEATUREIN.loc[:,['chr','start','end','sd_max','win_sim','sd_var','entropy','TEgini','sim_bubble','com_bubble','vntr_cov','vntr_gini_max']]
APG2=FEATUREIN.copy()
APG2['Vscore']=(1-FEATUREIN['win_sim'])+FEATUREIN['entropy']
APG2 = APG2.dropna()

APGcopy2=APG2.copy()
APGbef2=APGcopy2.iloc[:,[0,1,2]]
APGin2=APGcopy2[['sd_max','sd_var','sim_bubble','com_bubble','Vscore','TEgini','vntr_cov','vntr_gini_max']]
APGin2 = APGin2.dropna()

X = APGin2
import joblib
scaler = joblib.load(Modeldir+'/clustering_model/model/models/standard_scaler.pkl')  # 替换为你的pkl文件路径
X_scaled = scaler.transform(X)




def autoencoder(dims, act='relu', init='glorot_uniform'):
    """
    Fully connected auto-encoder model, symmetric.
    Arguments:
        dims: list of number of units in each layer of encoder. dims[0] is input dim, dims[-1] is units in hidden layer.
            The decoder is symmetric with encoder. So number of layers of the auto-encoder is 2*len(dims)-1
        act: activation, not applied to Input, Hidden and Output layers
    return:
        (ae_model, encoder_model), Model of autoencoder and model of encoder
    """
    n_stacks = len(dims) - 1
    input_img = Input(shape=(dims[0],), name='input')
    x = input_img
    # internal layers in encoder
    for i in range(n_stacks-1):
        x = Dense(dims[i + 1], activation=act, kernel_initializer=init, name='encoder_%d' % i)(x)
    # hidden layer
    encoded = Dense(dims[-1], kernel_initializer=init, name='encoder_%d' % (n_stacks - 1))(x)  # hidden layer, features are extracted from here
    x = encoded
    # internal layers in decoder
    for i in range(n_stacks-1, 0, -1):
        x = Dense(dims[i], activation=act, kernel_initializer=init, name='decoder_%d' % i)(x)
    # output
    x = Dense(dims[0], kernel_initializer=init, name='decoder_0')(x)
    decoded = x
    return Model(inputs=input_img, outputs=decoded, name='AE'), Model(inputs=input_img, outputs=encoded, name='encoder')

class ClusteringLayer(Layer):
    """
    Clustering layer converts input sample (feature) to soft label, i.e. a vector that represents the probability of the
    sample belonging to each cluster. The probability is calculated with student's t-distribution.

    # Example
    ```
        model.add(ClusteringLayer(n_clusters=10))
    ```
    # Arguments
        n_clusters: number of clusters.
        weights: list of Numpy array with shape `(n_clusters, n_features)` witch represents the initial cluster centers.
        alpha: degrees of freedom parameter in Student's t-distribution. Default to 1.0.
    # Input shape
        2D tensor with shape: `(n_samples, n_features)`.
    # Output shape
        2D tensor with shape: `(n_samples, n_clusters)`.
    """
    def __init__(self, n_clusters, weights=None, alpha=1.0, **kwargs):
        if 'input_shape' not in kwargs and 'input_dim' in kwargs:
            kwargs['input_shape'] = (kwargs.pop('input_dim'),)
        super(ClusteringLayer, self).__init__(**kwargs)
        self.n_clusters = n_clusters
        self.alpha = alpha
        self.initial_weights = weights
        self.input_spec = InputSpec(ndim=2)
    def build(self, input_shape):
        assert len(input_shape) == 2
        input_dim = input_shape[1]
        self.input_spec = InputSpec(dtype=K.floatx(), shape=(None, input_dim))
        self.clusters = self.add_weight((self.n_clusters, input_dim), initializer='glorot_uniform', name='clusters')
        if self.initial_weights is not None:
            self.set_weights(self.initial_weights)
            del self.initial_weights
        self.built = True
    def call(self, inputs, **kwargs):
        """ student t-distribution, as same as used in t-SNE algorithm.
         Measure the similarity between embedded point z_i and centroid µ_j.
                 q_ij = 1/(1+dist(x_i, µ_j)^2), then normalize it.
                 q_ij can be interpreted as the probability of assigning sample i to cluster j.
                 (i.e., a soft assignment)
        Arguments:
            inputs: the variable containing data, shape=(n_samples, n_features)
        Return:
            q: student's t-distribution, or soft labels for each sample. shape=(n_samples, n_clusters)
        """
        q = 1.0 / (1.0 + (tf.reduce_sum(tf.square(tf.expand_dims(inputs, axis=1) - self.clusters), axis=2) / self.alpha))
        q **= (self.alpha + 1.0) / 2.0
        q = tf.transpose(tf.transpose(q) / tf.reduce_sum(q, axis=1)) # Make sure each sample's 10 values add up to 1.
        return q
    def compute_output_shape(self, input_shape):
        assert input_shape and len(input_shape) == 2
        return input_shape[0], self.n_clusters
    def get_config(self):
        config = {'n_clusters': self.n_clusters}
        base_config = super(ClusteringLayer, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))


##run
x=X_scaled
from tensorflow.keras.optimizers import Adam
dims = [x.shape[-1], 128,128, 256,10]
init1 = VarianceScaling(scale=1. / 3., mode='fan_in',distribution='uniform')

autoencodernow, encodernow = autoencoder(dims, init=init1)
filename=Modeldir+"/clustering_model/model/models/dc.end.weights.h5"
autoencodernow.load_weights(filename)
pre=encodernow.predict(x)

from keras.models import Model
num_cluster=8
clustering_layer = ClusteringLayer(num_cluster, name='clustering')(encodernow.output)
model2 = Model(inputs=encodernow.input, outputs=[clustering_layer,autoencodernow.output])

filename = Modeldir+"/clustering_model/model/models/dc.fin.weights.h5"
model2.load_weights(filename)

q = model2.predict(x, verbose=0)
q=q[0]
y_pred = q.argmax(1)

complexnowdrow=APGcopy2.copy()
complexnowdrow['cluster']=y_pred
complexnowdrow_reset = complexnowdrow.reset_index(drop=True)

qx=q
n_clusters = qx.shape[1]  # 假设为 8
cluster_columns = [f'Cc{i}' for i in range(n_clusters)]

for i, col_name in enumerate(cluster_columns):
    complexnowdrow_reset[col_name] = qx[:, i]




########drawing....
out=complexnowdrow_reset.loc[:,['chr','start','end','cluster','Cc0','Cc2','Cc5','Cc3','Cc1','Cc4','Cc6','Cc7']]
replacement_dict = {1: "CC4", 2: "CC1", 0: "simple", 3: "CC3",
                   4: "CC5", 5: "CC2", 6: "CC6", 7: "CC7"}
out['cluster'] = out['cluster'].replace(replacement_dict)
out = out.rename(columns={
    'Cc0': 'simple',
    'Cc2': 'CC1',
    'Cc5': 'CC2',
    'Cc3':'CC3','Cc1':'CC4','Cc4':'CC5','Cc6':'CC6','Cc7':'CC7'
})




out.to_csv("final.anno",sep="\t",index=False)


out['POS']=(out['start']+out['end'])/2
cluster_columns = ['simple','CC1','CC2','CC3','CC4','CC5','CC6','CC7']
df = out[['POS'] + cluster_columns]  

colors=["#efefef",'#a87fb7','#d52918' ,'#2876bc',  '#a02735', '#684f7f','#da8e3d','#b96289']

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
from time import time
import numpy as np
import keras.backend as K
from tensorflow.keras.layers import Layer, InputSpec #from keras.engine.topology import Layer, InputSpec
from keras.layers import Dense, Input
from keras.models import Model
from keras.optimizers import SGD
from keras import callbacks
from keras.initializers import VarianceScaling
from sklearn.cluster import KMeans
import tensorflow as tf
import argparse
from scipy.interpolate import make_interp_spline
plt.figure(figsize=(35, 5))
for i, feature in enumerate(df.columns[1:]):
    x_new = np.linspace(df['POS'].min(), df['POS'].max(), 300)
    spline = make_interp_spline(df['POS'], df[feature], k=2)  # Cubic spline
    y_new = spline(x_new)
    y_new = np.clip(spline(x_new), 0, 1)
    plt.fill_between(x_new, y_new, color=colors[i], alpha=0.7, label=feature)  # 使用较低的透明度填充
    plt.plot(x_new, y_new, color=colors[i], linewidth=3)  # 绘制平滑曲线

plt.title('Genomic Region complexity annotation', fontsize=16)
plt.xlabel('POS', fontsize=14)
plt.ylabel('Probability', fontsize=14)
plt.xticks(np.arange(df['POS'].min(), df['POS'].max(), 100000), fontsize=12)

plt.yticks(fontsize=12)

plt.ylim(0, 1.1)
plt.legend(title='Cluster', fontsize=12)
plt.grid(visible=False)

plt.tight_layout()


plt.savefig('final_anno.complex.pdf', format='pdf', dpi=300)
plt.savefig('final_anno.complex.png', format='png', dpi=300)

#plt.show()