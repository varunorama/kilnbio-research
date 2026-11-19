"""Offline Circos-style views of explicitly mapped multi-omics measurements."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as DrawPath
from matplotlib.cm import ScalarMappable
from plot_design import CMAP, MISSING, INK, MUTED, PALETTE


def load_config(path):
    path=Path(path).resolve(); cfg=json.loads(path.read_text()); base=path.parent
    if cfg.get('schema_version')!=1 or cfg.get('coordinate_system')!='0-based-half-open':
        raise ValueError('Use schema_version 1 and coordinate_system 0-based-half-open')
    build=cfg.get('genome_build')
    if not isinstance(build,str) or not build.strip():raise ValueError('Declare a genome_build')
    files=[path]
    def read(name,required,limit):
        p=(base/name).resolve(); t=pd.read_csv(p,sep='\t',nrows=limit+1,dtype={'chrom':str,'chrom1':str,'chrom2':str})
        if len(t)>limit:raise ValueError(f'{p.name}: exceeds display limit {limit}; explicitly select or aggregate input first')
        if not required<=set(t):raise ValueError(f'{p.name}: requires columns {sorted(required)}')
        files.append(p); return t
    sizes=read(cfg['chrom_sizes'],{'chrom','length'},64)
    sizes['length']=pd.to_numeric(sizes.length,errors='raise')
    if sizes.empty or sizes.chrom.isna().any() or sizes.chrom.duplicated().any() or not np.isfinite(sizes.length).all() or (sizes.length<=0).any() or (sizes.length%1!=0).any():
        raise ValueError('Chromosome names must be unique and lengths positive integers')
    lengths=sizes.set_index('chrom').length.to_dict()
    def coordinates(t,suffix=''):
        chrom='chrom'+suffix;start='start'+suffix;end='end'+suffix
        for c in [start,end]:t[c]=pd.to_numeric(t[c],errors='raise')
        valid=t[chrom].isin(lengths) & np.isfinite(t[start]) & np.isfinite(t[end])
        valid &= (t[start]%1==0)&(t[end]%1==0)&(t[start]>=0)&(t[start]<t[end])&(t[end]<=t[chrom].map(lengths))
        if not valid.all():raise ValueError('Unknown chromosome or invalid/out-of-bounds interval; check genome build and coordinate convention')
    specs=cfg.get('tracks',[])
    if not isinstance(specs,list) or not 1<=len(specs)<=8:raise ValueError('Supply 1–8 measured tracks')
    names=[];tracks=[]
    for spec in specs:
        if spec.get('genome_build')!=build:raise ValueError('All tracks must explicitly match the declared genome build; no automatic liftover')
        name=spec.get('name');unit=spec.get('value_label')
        if not isinstance(name,str) or not name.strip() or name in names or not isinstance(unit,str) or not unit.strip():
            raise ValueError('Track names must be unique and each track needs a value_label')
        names.append(name)
        lo=float(spec['value_min']);hi=float(spec['value_max'])
        if not np.isfinite([lo,hi]).all() or lo>=hi:raise ValueError('Track limits must be finite and increasing')
        t=read(spec['file'],{'chrom','start','end','value'},50000);coordinates(t)
        raw=t.value.copy();t['value']=pd.to_numeric(raw,errors='raise')
        if np.isinf(t.value).any() or (t.value.dropna()<lo).any() or (t.value.dropna()>hi).any():
            raise ValueError('Measured values must be finite or missing and lie within declared track limits')
        t['track']=name;t['value_label']=unit
        tracks.append((spec,t))
    links=pd.DataFrame(columns=['chrom1','start1','end1','chrom2','start2','end2','relation'])
    spec=cfg.get('links')
    if spec is not None:
        if spec.get('genome_build')!=build:raise ValueError('Links must match the declared genome build')
        links=read(spec['file'],set(links.columns),200)
        coordinates(links,'1');coordinates(links,'2')
        if links.relation.isna().any() or links.relation.astype(str).str.strip().eq('').any():raise ValueError('Each link requires a relation label')
        if links.relation.nunique()>8:raise ValueError('At most eight link relation types can be displayed legibly')
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return cfg,sizes,tracks,links,hashes


def add_circos(report,config,allowed_layers=None):
    cfg,sizes,tracks,links,hashes=load_config(config)
    if allowed_layers is not None and not {s['name'] for s,_ in tracks}<=set(allowed_layers):
        raise ValueError('Circos track names must match layers in this multiomic-meta manifest')
    gap=min(.035,.8/len(sizes));span=2*np.pi-gap*len(sizes)
    starts={};widths={};offset=0
    for row in sizes.itertuples():
        width=span*row.length/sizes.length.sum();starts[row.chrom]=offset;widths[row.chrom]=width;offset+=width+gap
    lengths=sizes.set_index('chrom').length.to_dict()
    def angle(chrom,pos):return np.pi/2-(starts[chrom]+widths[chrom]*pos/lengths[chrom])
    def xy(theta,radius):return np.array([np.cos(theta),np.sin(theta)])*radius
    fig=plt.figure(figsize=(13,10),layout='constrained')
    grid=fig.add_gridspec(1,2,width_ratios=[3.7,1.55],wspace=.05)
    ax=fig.add_subplot(grid[0]);key=fig.add_subplot(grid[1]);key.axis('off')
    radii=np.linspace(.91,.47,len(tracks)) if len(tracks)>1 else np.array([.85])
    for i,row in enumerate(sizes.itertuples()):
        a=angle(row.chrom,0);z=angle(row.chrom,row.length)
        ax.add_patch(Wedge((0,0),1.02,np.degrees(z),np.degrees(a),width=.035,facecolor=PALETTE[i%len(PALETTE)],edgecolor='white',lw=.5))
        mid=(a+z)/2;pt=xy(mid,1.10)
        ax.text(*pt,row.chrom,ha='center',va='center',fontsize=9)
        for radius in radii:
            th=np.linspace(z,a,120);ax.plot(radius*np.cos(th),radius*np.sin(th),color='#E8EDF2',lw=9,zorder=0)
    values=[]
    key.text(0,1,'TRACKS / OUTSIDE → INSIDE',fontsize=10,fontweight='bold',color=INK,va='top')
    spacing=min(.092,.66/len(tracks))
    for i,((spec,t),radius) in enumerate(zip(tracks,radii)):
        theta=np.array([angle(r.chrom,(r.start+r.end)/2) for r in t.itertuples()])
        norm=Normalize(spec['value_min'],spec['value_max'])
        valid=t.value.notna().to_numpy()
        ax.scatter(radius*np.cos(theta[valid]),radius*np.sin(theta[valid]),c=t.loc[valid,'value'],norm=norm,cmap=CMAP,s=11,linewidths=0,zorder=3)
        ax.scatter(radius*np.cos(theta[~valid]),radius*np.sin(theta[~valid]),c=MISSING,marker='x',s=14,linewidths=.7,zorder=3)
        y=.93-i*spacing
        key.text(0,y,f'{i+1:02d}  {spec["name"]}',fontsize=11,fontweight='bold',va='top')
        key.text(0,y-.025,f'{spec["value_label"]} · {int(valid.sum()):,} measured · {int((~valid).sum()):,} missing',fontsize=8,color=MUTED,va='top')
        cax=key.inset_axes([0,y-.064,.84,.013])
        cb=fig.colorbar(ScalarMappable(norm=norm,cmap=CMAP),cax=cax,orientation='horizontal',ticks=[spec['value_min'],spec['value_max']])
        cb.ax.tick_params(labelsize=7,pad=1);cb.outline.set_visible(False)
        values.append(t.assign(display_radius=radius,display_angle_radians=theta))
    relation_colors={name:PALETTE[i%len(PALETTE)] for i,name in enumerate(sorted(links.relation.unique()))}
    link_radius=float(radii.min()-.07)
    for row in links.itertuples():
        p1=xy(angle(row.chrom1,(row.start1+row.end1)/2),link_radius)
        p2=xy(angle(row.chrom2,(row.start2+row.end2)/2),link_radius)
        path=DrawPath([p1,p1*.25,p2*.25,p2],[DrawPath.MOVETO,DrawPath.CURVE4,DrawPath.CURVE4,DrawPath.CURVE4])
        ax.add_patch(PathPatch(path,facecolor='none',edgecolor=relation_colors[row.relation],alpha=.55,lw=.9,zorder=1))
    y=.93-len(tracks)*spacing-.035
    key.text(0,y,f'DECLARED LINKS / n = {len(links):,}',fontsize=10,fontweight='bold')
    for j,(name,color) in enumerate(relation_colors.items()):
        key.plot([0,.08],[y-.035-j*.027]*2,color=color,lw=2)
        key.text(.11,y-.035-j*.027,name,fontsize=8,va='center')
    ax.text(0,0,f'{len(tracks)} assay tracks\n{len(links)} supplied links',ha='center',va='center',fontsize=11,color=MUTED,bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=8))
    key.set_xlim(0,1);key.set_ylim(0,1)
    ax.set_aspect('equal');ax.set_xlim(-1.22,1.22);ax.set_ylim(-1.22,1.22);ax.axis('off')
    fig.suptitle('Mapped multi-omics in genomic context\n'+cfg['genome_build']+' · independent measurement scales per track',x=.02,ha='left',fontsize=16,fontweight='bold')
    report.table('circos_links',links)
    report.table('circos_chromosomes',sizes)
    report.table('circos_inputs',pd.DataFrame([dict(input=k,sha256=v) for k,v in hashes.items()]))
    report.add('circos_tracks',fig,
        'Circos-style plot rendered by KilnBio/Matplotlib. Chromosome spans follow supplied lengths and order; '
        'dots mark interval midpoints. Rings follow configured assay order and retain independent declared value '
        'scales. Gray crosses denote supplied missing values; blank locations have no supplied row. Curves show '
        'only declared links and their relation labels, not interactions inferred from proximity. Overlapping dots '
        'may obscure each other; all input rows remain in TSV. No new test or target selection is performed.',pd.concat(values,ignore_index=True))
    return dict(genome_build=cfg['genome_build'],coordinate_system=cfg['coordinate_system'],
                tracks=[dict(name=s['name'],rows=len(t),value_label=s['value_label'],value_min=s['value_min'],value_max=s['value_max']) for s,t in tracks],
                links=len(links),inputs_sha256=hashes,renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                citation='Krzywinski et al. 2009. Genome Research 19:1639–1645. doi:10.1101/gr.092759.109; circular visualization concept, not the original Circos executable.')


def run(config,out):
    from visual_reports import Report
    out=Path(out).resolve()
    if out.exists():raise ValueError('Use a new output directory')
    load_config(config)
    out.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='.circos-',dir=out.parent))
    try:
        report=Report(stage,'KilnBio circular multi-omics report',[
            'Measurements and links are supplied explicitly. All tracks must declare the same genome build.',
            'No missing assay is fabricated. Protein/RNA values require reviewed mapping to genomic intervals.'])
        receipt=add_circos(report,config);report.finish(receipt)
        if out.exists():raise ValueError('Output directory appeared during rendering')
        stage.rename(out)
    finally:
        if stage.exists():shutil.rmtree(stage)
    return out/'report.html'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('config');p.add_argument('out')
    a=p.parse_args()
    try:print(run(a.config,a.out))
    except (ValueError,OSError,KeyError,TypeError) as e:p.exit(2,str(e)+'\n')


if __name__=='__main__':main()
