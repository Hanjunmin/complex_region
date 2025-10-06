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
parser.add_argument('--e', help='epochs')
parser.add_argument('--p', help='Dimensionality of dimensionality reduction')
args = parser.parse_args()
epo=int(args.e)
dimen=int(args.p)



initdata2 = pd.read_csv("./COMPLEXFEATURE/05_split/allfeaturen",sep=",",header=0)  
# initdata2.columns=['chr','start','end','SDsim','vntr','Vscore','entropy','Vscoren','trimTEsd','TEmax','trimsdout']
initdata2.columns=['chr','start','end','entropy', 'win_sim','TEvar','TEgini','sd_max' ,   'sd_var']
 

###添加VNTR
VNTR= pd.read_csv("VNTR.sd.cov",header=0,sep="\t")  ##allfeaturen是后来改了vscore的结果
initdata2=pd.merge(initdata2, VNTR, how='left', on=['chr','start','end'])
initdata2.loc[initdata2['vntr_cov']>0,"vntr_cov"]=1
# initdata2.loc[initdata2['vntr_sd_max']>0.1,"vntr_sd_max"]=1

###添加sd
# SD= pd.read_csv("AL_SD.sim_sd",header=None,sep="\t")  ##allfeaturen是后来改了vscore的结果
# SD.columns=['chr','start','end','win','SDsimi','SD_sd']
# SD['start']=SD['start']-1
# SD['SD_sd']=SD['SD_sd']*SD['SDsimi']

#initdata2=pd.merge(initdata2, SD, how='left', on=['chr','start','end'])



initdata2['sim_bubble']=0.01
initdata2['com_bubble']=0.01
refinebubble = pd.read_csv("./COMPLEXFEATURE/09_bubble_only/02_INTE/end.txt",header=None,sep="\t")
refinebubble.columns=['chr','start','end','chrregion','startregion','endregion','poly','polyratio','simple','pathratio','lenration','pos','polyrefine','H','V4new','V5new','norm_bub','bub']
refinebubble['sim_bub']=refinebubble['V5new']+refinebubble['H']
refinebubble['com_bub']=refinebubble['V4new']+refinebubble['H']
refinebubble['start']=refinebubble['start']-1
initdatac=pd.merge(initdata2, refinebubble, how='left', on=['chr','start','end'])
initdatac.loc[initdatac['sim_bub'].notna(), 'sim_bubble'] = initdatac['sim_bub']
initdatac.loc[initdatac['com_bub'].notna(), 'com_bubble'] = initdatac['com_bub']




# refinebubble = pd.read_csv("./COMPLEXFEATURE/09_bubble_only/02_INTE/end.txt",header=None,sep="\t")
# if bubty==1: ##51014 490351
#     refinebubble.columns=['chr','start','end','chrregion','startregion','endregion','poly','polyratio','simple','pathratio','lenration','pos','polyrefine','H','V4new','V5new','norm_bub','bub']
# else:   ##51003  490442
#     refinebubble.columns=['chr','start','end','chrregion','startregion','endregion','poly','polyratio','simple','pathratio','lenration','pos','polyrefine','H','V4new','V5new','bub','norm_bub']



initdatac=initdatac.loc[:,['chr','start','end','sd_max','win_sim','sd_var','entropy','TEgini','sim_bubble','com_bubble','vntr_cov','vntr_gini_max']]
# APG=initdata[initdata['chr']=="chr1"]
APG2=initdatac.copy()
APG2['Vscore']=(1-initdatac['win_sim'])+initdatac['entropy']
APG2 = APG2.dropna()

# APG2.loc[APG2['vntr'] < 0.05, 'vntr'] = 0
# APG2.loc[APG2['simple'] <30 , 'insertion'] = 0
# APG2.loc[APG2['simple'] <30 , 'super'] = 0
# APG2.loc[APG2['simple'] <30 , 'simple'] = 0
# APG2.loc[APG2['Vscore'] <0.14 , 'Vscore'] = 0
APGcopy2=APG2.copy()
APGbef2=APGcopy2.iloc[:,[0,1,2]]
APGin2=APGcopy2[['sd_max','sd_var','sim_bubble','com_bubble','Vscore','TEgini','vntr_cov','vntr_gini_max']]
APGin2 = APGin2.dropna()

X = APGin2
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)  # 标准化处理




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
dims = [x.shape[-1], 128,128, 256, dimen]
init1 = VarianceScaling(scale=1. / 3., mode='fan_in',distribution='uniform')
pretrain_optimizer = SGD(learning_rate=0.05 ,momentum=0.9)
#pretrain_optimizer = Adam(learning_rate=0.01)
pretrain_epochs =  epo ###300
batch_size = 3000  ###300


##run
autoencodernow, encodernow = autoencoder(dims, init=init1)
# feature=4
# for i in range(10):
autoencodernow.compile(optimizer=pretrain_optimizer, loss='mse') 
autoencodernow.fit(x, x, batch_size=batch_size, epochs=pretrain_epochs) #, callbacks=cb)
# filename="iter"+str(i)+"_"+str(feature)+".weights.h5"
filename="lianxu"+str(epo)+"_"+str(dimen)+"_dc.end.weights.h5"
autoencodernow.save_weights(filename)


##run
autoencodernow, encodernow = autoencoder(dims, init=init1)
filename="lianxu"+str(epo)+"_"+str(dimen)+"_dc.end.weights.h5"
autoencodernow.load_weights(filename)
pre=encodernow.predict(x)
# delete_indices = complexnow[complexnow['simple'] >20000].index
# clusterxinenout=np.delete(pre, delete_indices,axis=0)
# clusterxin=np.delete(x, delete_indices,axis=0)
##run
# from tensorflow.keras.optimizers import Adam
# num_cluster=5
# clustering_layer = ClusteringLayer(num_cluster, name='clustering')(encodernow.output)
# model = Model(inputs=encodernow.input, outputs=clustering_layer)
# for layer in model.layers[:-1]:  # 假设聚类层是最后一层
#     layer.trainable = False

# model.layers[-1].trainable = True  # 聚类层
# model.compile(optimizer=Adam(learning_rate=0.001), loss='kld')


# ##run
# complexnow=APGcopy2.copy()

# delete_indices = complexnow[complexnow['simple'] >30000].index
# clusterxinenout=np.delete(pre, delete_indices,axis=0)
# clusterxin=np.delete(x, delete_indices,axis=0)

num_cluster=8
kmeans = KMeans(n_clusters=num_cluster, n_init=20)
y_pred = kmeans.fit_predict(pre)
# complexnowdrow=complexnow[complexnow['simple'] <=30000]
# complexnowdrow['cluster']=y_pred
# complexnowdrow_reset = complexnowdrow.reset_index(drop=True)
APGcopy2['cluster']=y_pred
filename="lianxu"+str(epo)+"_"+str(dimen)+"dc._csv"

APGcopy2.to_csv(filename)


# import pandas as pd
# import numpy as np
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler
# import seaborn as sns
# import matplotlib.pyplot as plt
# from sklearn.mixture import GaussianMixture
# plt.figure(figsize=(14, 8))
# for i, feature in enumerate(['Vscore', 'SDsim', 'simple', 'super', 'insertion','vntr']):
#     plt.subplot(2, 3, i+1)
#     sns.boxplot(x='cluster', y=feature, data=complexnowdrow)
#     plt.title(f'Distribution of {feature} by cluster')


# plt.tight_layout()
# plt.show()