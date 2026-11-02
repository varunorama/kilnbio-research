"""Plot completed control experiments; output stays beside the local results."""
import json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main(base):
    s=json.loads((base/'summary.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8.7),layout='constrained')
    fig.suptitle('KilnBio: biological controls across assay entry points',fontsize=17,fontweight='bold')
    ax=axes[0,0];records=s['airway']+s['kang'];labels=['CRISPLD2','IFIT1','ISG15','MX1']
    for engine,color,delta in [('deseq2','#286ca8',-.1),('edger','#d97935',.1)]:
        ys=[next(r['log2_fold_change'] for r in records if r['engine']==engine and r['symbol']==g) for g in labels]
        ax.scatter(np.arange(4)+delta,ys,color=color,s=50,label={'deseq2':'DESeq2','edger':'edgeR'}[engine])
    ax.axhline(0,color='.7',linewidth=.7);ax.set_xticks(range(4),['CRISPLD2\n4 paired donors','IFIT1\n8 paired donors','ISG15\n8 paired donors','MX1\n8 paired donors'])
    ax.set_ylabel('Model log2 fold change');ax.set_title('A  Published stimulus-response controls',loc='left',fontweight='bold')
    donor_note='Every donor increased for each marker' if all(r['donors_up']==r['donors'] for r in records) else 'Mixed individual-donor directions; see summary'
    ax.legend(frameon=False);ax.text(.02,.98,'Bulk RNA: dexamethasone | 10x monocytes: IFN-beta\n'+donor_note,transform=ax.transAxes,va='top',fontsize=9)
    ax.set_ylim(-.4,11)
    ax=axes[0,1];r=s['atac_dnase'];v=[100*x['direction_agreement'] for x in r]
    ax.bar([x['cell'] for x in r],v,color=['#399a88','#286ca8'],width=.55)
    for i,x in enumerate(r):ax.text(i,v[i]+2,f'{v[i]:.1f}%\n{x["tested_tiles"]:,} tested tiles',ha='center',fontsize=9)
    ax.set_ylim(0,122);ax.set_ylabel('ATAC agrees with cell-specific DNase (%)')
    ax.set_title('B  Accessibility supported by another assay',loc='left',fontweight='bold')
    ax.text(.02,.97,'Fixed 10 kb autosomal tiles; 2 cultures per cell line\nRegional agreement is not a donor-level accuracy',transform=ax.transAxes,va='top',fontsize=9)
    ax=axes[1,0];r=s['colon']['paired_common_coverage'];means=[100*x['mean_cancer_minus_normal'] for x in r]
    medians=[100*x['median_cancer_minus_normal'] for x in r]
    ax.scatter(range(1,4),means,label='Mean',s=60,color='#286ca8');ax.scatter(range(1,4),medians,label='Median',marker='s',s=45,color='#d97935')
    ax.axhline(0,color='.6',linewidth=.8);ax.set_xticks([1,2,3],['Patient 1','Patient 2','Patient 3']);ax.set_ylim(-20,5)
    ax.set_ylabel('Tumor minus normal methylation (pp)');ax.legend(frameon=False,loc='lower left')
    ax.set_title('C  Paired WGBS methylation shifts',loc='left',fontweight='bold')
    ax.text(.02,.97,f'{r[0]["CpGs"]:,} chr21 CpGs with depth >=5 in all six samples\nMeans and medians shown separately for each patient',transform=ax.transAxes,va='top',fontsize=9)
    ax=axes[1,1]
    for i,cell in enumerate(['GM12878','K562']):
        y=[next(r['DHS_minus_absent_median_log2_enrichment'] for r in s['chip_input_adjusted_exploratory'] if r['cell']==cell and r['assay']==mark) for mark in ['H3K27me3','H3K36me3']]
        ax.bar(np.arange(2)+(i-.5)*.32,y,width=.30,label=cell,color=['#286ca8','#399a88'][i])
    ax.set_xticks([0,1],['H3K27me3','H3K36me3']);ax.axhline(0,color='.6',linewidth=.8);ax.set_ylim(-.9,.9)
    ax.set_ylabel('DHS minus absent-DHS log2 IP/input enrichment')
    ax.set_title('D  Histone context against independent DNase',loc='left',fontweight='bold')
    ax.legend(frameon=False,loc='lower right');ax.text(.02,.97,'Exploratory, pooled cultures; 10 kb resolution\nAssess each mark and cell line separately',transform=ax.transAxes,va='top',fontsize=9)
    fig.savefig(base/'biological-controls.png',dpi=180);fig.savefig(base/'biological-controls.svg');plt.close(fig)

if __name__=='__main__':main(Path(sys.argv[1]).resolve())
