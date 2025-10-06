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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import  umap
from sklearn.manifold import TSNE
# parser = argparse.ArgumentParser(description='RUNMODEL')
# parser.add_argument('--s', help='sample number')
# parser.add_argument('--nb', help='neighbors')
# parser.add_argument('--d', help='min dist')
# args = parser.parse_args()
# samnum=int(args.s)
# nb=int(args.nb)
# dis=float(args.d)


###chr1:12054000-12954000
#cluster_data = pd.read_csv("./COMPLEXFEATURE/09_bubble_only/03_model/APGcopy2.csv",header=0)  
cluster_data = pd.read_csv("./COMPLEXFEATURE/09_bubble_only/03_model/lianxu150_10_final.csv",header=0) 


filein="./COMPLEXFEATURE/09_bubble_only/03_model/lianxu150_10_final.csv"
APGcopy2 = pd.read_csv(filein,header=0)  
column_mapping = {
    'sd_max': 'inter_sim',
    'sd_var': 'inter_var',
    'TEgini': 'TE_var',
    'sim_bubble': 'global_bubble',
    'com_bubble': 'local_bubble',
    'vntr_cov': 'intra_sim',
    'vntr_gini_max': 'intra_var',
    'Vscore': 'seq_div'
}


df_renamed = APGcopy2.rename(columns=column_mapping)
cluster_mapping = {
    0: 'simple',
    1: 'CC4', 
    2: 'CC1',
    3: 'CC3',
    4: 'CC5',
    5: 'CC2',
    6: 'CC6',
    7: 'CC7',
}
df_renamed['cluster'] = df_renamed['cluster'].replace(cluster_mapping)
cluster_data=df_renamed

cluster_datacom=cluster_data[cluster_data['cluster']!="simple"]
cluster_datacom=cluster_datacom.sample(n=10000, random_state=42)

cluster_datasim=cluster_data[cluster_data['cluster']=="simple"]
cluster_datasim= cluster_datasim.sample(n=10000, random_state=42)  
cluster_dataother=cluster_data[(cluster_data['chr']=="chr1") & (cluster_data['start']>=12054000) & (cluster_data['end']<=12954000)]

cluster_datacom=pd.concat([cluster_datacom,cluster_datasim,cluster_dataother])
cluster_datacom = cluster_datacom.drop_duplicates()
cluster_datacom = cluster_datacom.reset_index(drop=True)  



selected_features=cluster_datacom[['inter_sim','inter_var','global_bubble','local_bubble','seq_div','intra_sim','intra_var','TE_var' ]]
X = selected_features.values
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
sample_indices=np.array(cluster_datacom.index)
X_sampled = X_scaled[sample_indices]
cluster_labels = cluster_datacom['cluster']
X_sampled = X_scaled[sample_indices]
labels_sampled = cluster_labels[sample_indices]
subdata = selected_features.iloc[sample_indices,:]
from sklearn.preprocessing import LabelEncoder
le = LabelEncoder()
labels_sampled_encoded = le.fit_transform(labels_sampled)

unique_clusters = np.unique(labels_sampled_encoded)
colors = plt.cm.get_cmap('tab20', len(unique_clusters))


import numpy as np
import matplotlib.pyplot as plt
features_to_visualize = ['inter_sim','inter_var','global_bubble','local_bubble','seq_div','intra_sim','intra_var','TE_var' ]
# colormaps = ['viridis', 'plasma', 'magma', 'cividis', 
#              'inferno', 'coolwarm', 'Spectral', 'turbo']

n_features = len(features_to_visualize)
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


####umap
umap_model = umap.UMAP(
    n_components=2,          
    n_neighbors=300,          
    min_dist=1,           
    random_state=42          
)
X_umap = umap_model.fit_transform(X_sampled)



plt.close('all')
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))  # Width=4*4=16, Height=4*2=8
axes = axes.flatten()  # Flatten to easily iterate over 8 subplots
for i, (ax, feature_name) in enumerate(zip(axes, features_to_visualize)):
    if i >= 8:  # Break if we have more features than subplots
        break
    feature_idx = features_to_visualize.index(feature_name)
    feature_values = subdata.iloc[:, feature_idx]
    current_cmap = 'viridis'  # Good-looking colormap
    scatter = ax.scatter(
        X_umap[:, 0],
        X_umap[:, 1],
        c=feature_values,
        cmap=current_cmap,
        s=15,
        alpha=0.6,
        edgecolors='w',
        linewidths=0.3,
        vmin=feature_values.min(),
        vmax=feature_values.max()
    )
    ax.set_title(f"{feature_name}", pad=10, fontsize=10)
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.grid(True, linestyle='--', alpha=0.3)
    
    # Add colorbar to the right of each subplot (without label)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(scatter, cax=cax, ax=ax)
    cbar.ax.tick_params(labelsize=8)

for i in range(len(features_to_visualize), len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.savefig("supp_feature.png", format='png', bbox_inches='tight')
plt.savefig("supp_feature.pdf", format='pdf', bbox_inches='tight', dpi=300)

plt.show()








colorsdict = {
    'simple':'#A5A5A5',
    'CC4':'#a02735', 
    'CC1':'#a87fb7',
    'CC3':'#2876bc',
    'CC5':'#684f7f',
    'CC2':'#d52918',
    'CC6':'#da8e3d',
    'CC7':'#b96289',
}

import matplotlib
# Use a stable backend to avoid Tkinter issues
matplotlib.use('qt5agg')  # or 'agg' for non-interactive, or 'wxagg' if preferred
import matplotlib.pyplot as plt
import numpy as np

# Close any existing figures to prevent Tkinter conflicts
plt.close('all')

# Create a single figure
plt.figure(figsize=(12, 8))

# Plot scatter points for each cluster
for cluster_encoded, cluster_str in zip(unique_clusters, np.unique(labels_sampled)):
    mask = (labels_sampled_encoded == cluster_encoded)
    
    # Plot sample_indices1 (chr1 region): larger, more opaque, square markers
    mask1 = mask & np.isin(sample_indices, sample_indices1)
    if np.any(mask1):
        plt.scatter(
            X_umap[mask1, 0], X_umap[mask1, 1],
            color=colorsdict[cluster_str],
            marker='s',  # Square
            s=30,        # Larger points
            alpha=0.9,   # More opaque
            edgecolors='w',  # White edges
            linewidths=0.5
        )
    
    # Plot sample_indices2 (random sampling): smaller, more transparent, circular markers
    mask2 = mask
    if np.any(mask2):
        plt.scatter(
            X_umap[mask2, 0], X_umap[mask2, 1],
            color=colorsdict[cluster_str],
            marker='o',  # Circular
            s=10,        # Smaller points
            alpha=0.2,   # More transparent
            edgecolors='none'  # No edges
            # No label to avoid influencing legend
        )
    
    # Create a hidden scatter plot for solid legend markers
    plt.scatter(
        [], [],  # Empty data to not display on plot
        color=colorsdict[cluster_str],
        marker='o',  # Circular
        s=10,        # Base size (will be scaled in legend)
        alpha=1.0,   # Fully solid for legend
        edgecolors='none',
        label=f'{cluster_str}'  # Label for legend
    )

# Customize legend in the right upper corner with larger, solid markers
plt.legend(
    bbox_to_anchor=(1.05, 1),  # Right upper corner
    loc='upper left',
    fontsize=9,
    markerscale=2.0,  # Make legend markers ~2x larger than s=10
    handlelength=2,   # Length of legend markers
    handletextpad=0.5,  # Spacing between marker and text
    frameon=True,     # Add a frame around the legend
    labelcolor='black'  # Text color for readability
)

# Add axis labels and grid
plt.xlabel("UMAP 1", fontsize=14)
plt.ylabel("UMAP 2", fontsize=14)
plt.grid(alpha=0.1)

# Adjust layout to accommodate legend
plt.tight_layout()

# Save the plot
plt.savefig("1_24Mumap.pdf", format='pdf', bbox_inches='tight', dpi=300)
plt.savefig("1_24Mumap.png", format='png', bbox_inches='tight')

# Show the plot
plt.show()



