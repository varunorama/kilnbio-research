"""Create fictional RNA/Olink/methylation inputs and a single-command study recipe."""
import argparse
import json
from pathlib import Path
from make_synthetic import make


def create(out):
    out=Path(out).resolve();make(out)
    cfg=json.loads((out/'config.json').read_text())
    cfg['donors']='{{study}}/'+cfg['donors']
    for layer in cfg['layers']:
        layer['matrix']=('{{run}}/' if layer['name']=='olink' else '{{study}}/')+layer['matrix']
        layer['specimens']='{{study}}/'+layer['specimens']
        if layer['feature_map']!='identity':layer['feature_map']='{{study}}/'+layer['feature_map']
    (out/'layers-template.json').write_text(json.dumps(cfg,indent=2)+'\n')
    study={'schema_version':1,'research_use':True,'allow_remote_scientific_services':False,
        'inputs':['donors.tsv','specimens.tsv','olink_samples.tsv','rna.tsv','methylation.tsv','methylation_map.tsv','olink_map.tsv','olink_export.tsv'],
        'configs':{'layers':'layers-template.json'},'steps':[
        {'id':'olink','command':'olink-import','args':['{{study}}/olink_export.tsv','--samples','{{study}}/olink_samples.tsv','--out','{{run}}/olink'],
         'depends_on':[],'outputs':['olink/eligible_npx.tsv','olink/measurement_qc.tsv']},
        {'id':'layers','command':'multiomic-layers','args':['{{run}}/configs/layers.json','{{run}}/layers'],
         'depends_on':['olink'],'outputs':['layers/manifest.tsv','layers/layers_design.tsv','layers/rna.tsv','layers/olink.tsv','layers/methylation.tsv']},
        {'id':'joint','command':'multiomic-meta','args':['{{run}}/layers/manifest.tsv','--min-layers','2','--out','{{run}}/joint'],
         'depends_on':['layers'],'outputs':['joint/genes.tsv','joint/nominations.tsv','joint/summary.json']}]}
    (out/'study.json').write_text(json.dumps(study,indent=2)+'\n')
    print('Synthetic study created. Run: kilnbio study '+str(out/'study.json')+' NEW_RUN')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output');a=p.parse_args();create(a.output)
