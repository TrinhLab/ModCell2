import pandas as pd
import seaborn as sns
import matplotlib.patches as mpatches

figsize = (5,7)

df = pd.read_csv('fluxes.csv')
df3 = df.set_index('id')

# Categorical coloring
lut = dict(zip(df['subsystem'].unique(),sns.color_palette("husl", len(df['subsystem'].unique()))))
id_to_col = dict(zip(df3.index, df3['subsystem'].map(lut)))
row_colors = [id_to_col[id] for id in list(df3.index)]
df3.drop('subsystem','columns',inplace=True)

g = sns.clustermap(df3, cmap="RdBu",row_colors=row_colors, center=0, vmin=df3.values.min(), vmax=df3.values.max(), linewidths=.15, figsize=figsize,cbar_kws={"label": "Scaled flux"})

# categorical row legend
legend_TN = [mpatches.Patch(color=c, label=l) for l, c in lut.items()]
l2=g.ax_heatmap.legend(loc=(0,20),bbox_to_anchor=(1.01,1.01),handles=legend_TN,frameon=True)
l2.set_title(title='Subsystem',prop={'size':10})

#axis labels
ax = g.ax_heatmap
ax.set_xlabel('Production network')
ax.set_ylabel('Reaction')

# save
g.savefig('flux_cm.svg', format='svg')


## Turnover
figsize = (5,17)
df = pd.read_csv('turnover.csv')
df3 = df.set_index('id')

g = sns.clustermap(df3, cmap="RdBu", center=0, vmin=df3.values.min(), vmax=df3.values.max(), linewidths=.15, figsize=figsize,cbar_kws={"label": "Scaled turnover rate"})

ax = g.ax_heatmap

ax.set_xlabel('Production network')
ax.set_ylabel('Metabolite')

# save
g.savefig('turnover_cm.svg', format='svg')

## Turnover precursors and currency
figsize2 = (5,5)
df = pd.read_csv('turnover_prec.csv')
df3 = df.set_index('id')

g = sns.clustermap(df3, cmap="RdBu", center=0, vmin=df3.values.min(), vmax=df3.values.max(), linewidths=.15, figsize=figsize2,cbar_kws={"label": "Scaled turnover rate"})

ax = g.ax_heatmap

ax.set_xlabel('Production network')
ax.set_ylabel('Metabolite')

# save
g.savefig('turnover_prec_cm.svg', format='svg')
