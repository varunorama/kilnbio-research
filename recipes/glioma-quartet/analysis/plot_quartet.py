"""Show paired technical and biological errors, with the data kept outside Git."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
p=Path(sys.argv[1])
d=pd.read_csv(p/'held_out_labs.tsv',sep='\t')
fig,axes=plt.subplots(1,3,figsize=(10.6,3.9),layout='constrained')
for ax,assay in zip(axes,['RNA','protein','metabolite']):
    a=d[d.assay==assay]; before=a[a.method=='published'].set_index('lab')
    for name,color,marker,label in [('frozen_control_PC','#216a70','o','Frozen control projection'),('concurrent_D6','#ac5c24','D','Concurrent D6 reference')]:
        after=a[a.method==name].set_index('lab')
        delta=100*(after[['technical_RMSE','biological_contrast_RMSE']]/before[['technical_RMSE','biological_contrast_RMSE']]-1)
        ax.scatter(delta.technical_RMSE,delta.biological_contrast_RMSE,c=color,marker=marker,label=label,
                   alpha=.8,edgecolors='white',linewidths=.45,s=34,zorder=3)
    ax.axvline(0,color='0.55',lw=.8);ax.axhline(0,color='0.55',lw=.8)
    ax.set_title(f'{assay.capitalize() if assay!="RNA" else assay} · {len(before)} labs')
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlabel('Same-material error change (%)')
axes[0].set_ylabel('Biological contrast error change (%)')
axes[0].legend(fontsize=7,loc='best',frameon=False)
fig.suptitle('A reduction in technical error can accompany increased contrast error')
for ext in ['png','svg']:fig.savefig(p/('hidden-lab-tradeoff.'+ext),dpi=180)
plt.close(fig)
