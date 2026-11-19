"""Create fictional pathway evidence and a complete, runnable study; no downloads."""
import argparse
import json
from pathlib import Path
import csv


def create(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    sets={'SYNTHETIC_RESPONSE':['G1','G2','G3','G4'],'SYNTHETIC_OVERLAP':['G3','G4','G5'],
          'SYNTHETIC_DISCORDANT':['G6','G7'],'SYNTHETIC_NULL':['G8','G9'],'SYNTHETIC_UNMEASURED':['G10']}
    (out/'sets.gmt').write_text(''.join(n+'\tfictional teaching set\t'+'\t'.join(g)+'\n' for n,g in sets.items()))
    def write(name,headers,rows):
        with (out/name).open('w',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(headers);w.writerows(rows)
    layers=[]
    for i,layer in enumerate(['RNA','protein','chromatin']):
        write(layer+'_background.tsv',['gene_id'],[[f'G{g}'] for g in range(1,10)])
        write(layer+'.tsv',['pathway','p','score','contributing_genes'],[
            ['SYNTHETIC_RESPONSE',.0001*(i+1),2+i/4,'G1;G3' if i<2 else 'G4'],
            ['SYNTHETIC_OVERLAP',.03*(i+1),1.1,'G3;G5'],
            ['SYNTHETIC_DISCORDANT',.0002,(-1 if i==1 else 1)*2,'G6' if i==0 else 'G7'],
            ['SYNTHETIC_NULL',.8,1,''],['SYNTHETIC_UNMEASURED','NA','NA','']])
        layers.append(dict(layer=layer,results=layer+'.tsv',background=layer+'_background.tsv',method='Fictional valid-p teaching fixture, not a biological analysis',score_kind='signed',complete_family_reviewed=True))
    cfg=dict(gmt='sets.gmt',collection='Fictional teaching collection v1',species='human',gene_id_type='synthetic',contrast='fictional treated vs control',annotation_reviewed=True,min_assays=2,alpha=.05,joint_test=True,layers=layers)
    (out/'compare.json').write_text(json.dumps(cfg,indent=2)+'\n')
    # Study resolves input paths before dispatch; this also makes it MCP runnable.
    template=dict(cfg,gmt='{{study}}/sets.gmt',layers=[dict(x,results='{{study}}/'+x['results'],background='{{study}}/'+x['background']) for x in layers])
    (out/'compare.study.json').write_text(json.dumps(template,indent=2)+'\n')
    inputs=['sets.gmt',*[x[k] for x in layers for k in ['results','background']]]
    study=dict(schema_version=1,research_use=True,inputs=inputs,configs={'pathways':'compare.study.json'},steps=[dict(id='pathways',command='pathway-compare',args=['{{run}}/configs/pathways.json','{{run}}/pathways'],depends_on=[],outputs=['pathways/pathways.tsv','pathways/assay_evidence.tsv','pathways/report.html','pathways/visualizations.pdf','pathways/visualization.json'])])
    (out/'study.json').write_text(json.dumps(study,indent=2)+'\n')
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('output');args=parser.parse_args();print(create(args.output))
