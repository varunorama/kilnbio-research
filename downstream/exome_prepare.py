"""Public nf-core/sarek 3.12.0 adapter, scoped to autosomal paired-end DNA.

Public basis and the supported subset are documented in docs/EXOME_AND_CH.md.
"""
import csv
import re
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


def validate_params(assay, params):
    require(params.get('wes') is True, 'Set wes:true for an exome or capture panel')
    require(params.get('step', 'mapping') == 'mapping', 'This adapter starts from paired FASTQs')
    require(params.get('aligner') in ('bwa-mem', 'bwa-mem2', 'dragmap'), 'Choose an open DNA aligner explicitly')
    require(params.get('intervals') and params.get('fasta_fai'), 'Supply target BED and matched fasta_fai')
    sizes = {}
    with Path(params['fasta_fai']).open() as handle:
        for line in handle:
            row = line.split('\t')
            require(len(row) >= 2 and row[0] not in sizes, 'Invalid or duplicate FASTA index contig')
            sizes[row[0]] = int(row[1])
    intervals = 0
    with Path(params['intervals']).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                continue
            row = line.rstrip().split('\t')
            require(len(row) >= 3, 'Targets must be a tab-separated BED')
            chrom, start, end = row[0], int(row[1]), int(row[2])
            require(re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', chrom),
                    'Initial exome/CH route requires autosomal targets; sex-chromosome ploidy needs a dedicated protocol')
            require(chrom in sizes and 0 <= start < end <= sizes[chrom], 'Target outside reference or contig names disagree')
            intervals += 1
    require(intervals > 0, 'Target BED is empty')
    tools = params.get('tools', '').split(',')
    supported = {'haplotypecaller', 'deepvariant', 'strelka', 'manta', 'vep'} if assay == 'exome' else {'mutect2', 'vep'}
    require(bool(tools) and set(tools) <= supported and len(tools) == len(set(tools)), 'Choose supported tools explicitly')
    require(('mutect2' in tools) if assay == 'ch-somatic' else bool(set(tools) & {'haplotypecaller','deepvariant','strelka'}),
            'A small-variant caller must be selected')
    skipped = set(params.get('skip_tools', '').split(',')) - {''}
    require(skipped <= {'baserecalibrator','haplotypecaller_filter'}, 'Only explicit BQSR or HaplotypeCaller post-filter skipping is supported')
    require('haplotypecaller_filter' not in skipped or 'haplotypecaller' in tools,
            'haplotypecaller_filter skip requires that caller; resulting VCFs are unfiltered')
    if 'baserecalibrator' not in skipped:
        require(all(params.get(k) for k in ('dbsnp','dbsnp_tbi','known_indels','known_indels_tbi')),
                'BQSR requires local known sites and indexes, or an explicit skip_tools:baserecalibrator choice')
    if assay == 'ch-somatic':
        require(all(params.get(k) for k in ('pon','pon_tbi','germline_resource','germline_resource_tbi',
                                           'contamination_resource','contamination_resource_tbi')),
                'CH somatic calling requires a matched panel of normals, germline resource, common-SNP contamination resource and indexes')
    if 'vep' in tools:
        require(params.get('vep_cache') and '://' not in params['vep_cache'] and Path(params['vep_cache']).is_dir(),
                'VEP needs a local versioned cache')
        require(all(params.get(k) for k in ('vep_cache_version','vep_genome','vep_species','vep_custom_args')),
                'Declare VEP cache version, assembly, species and annotation arguments')
    for name in ('vep_spliceai','vep_dbnsfp','vep_mastermind','vep_condel','vep_phenotypes'):
        require(not params.get(name), f'{name} is outside the verified public-resource preset')
    params.setdefault('save_mapped', True)
    params.setdefault('save_output_as_bam', True)
    params.setdefault('normalize_vcfs', True)


def validate_samples(assay, sheet, params, local_file):
    expected = ['patient','sex','status','sample','lane','fastq_1','fastq_2']
    with sheet.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == expected, 'Sarek samplesheet requires ' + ','.join(expected))
        rows = list(reader)
    require(rows, 'Empty samplesheet')
    files, lanes, metadata = set(), set(), {}
    for row in rows:
        require(None not in row and all(v and v == v.strip() for v in row.values()), 'Missing or padded samplesheet values')
        for field in ('patient','sample','lane'):
            require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', row[field]), 'Unsafe DNA identifier')
        require(row['sex'] in ('XX','XY','NA'), 'sex must be XX, XY or NA')
        require(row['status'] == ('0' if assay == 'exome' else '1'),
                'Use status 0 for germline; CH blood uses status 1 solely for the somatic caller interface')
        signature = tuple(row[k] for k in ('patient','sex','status'))
        require(row['sample'] not in metadata or metadata[row['sample']] == signature, 'Conflicting sample metadata across lanes')
        metadata[row['sample']] = signature
        lane = (row['patient'], row['sample'], row['lane'])
        require(lane not in lanes, 'Duplicate lane for a sample')
        lanes.add(lane)
        for field in ('fastq_1','fastq_2'):
            path = local_file(row[field], sheet.parent)
            require(path.endswith(('.fq.gz','.fastq.gz')) and path not in files, 'Duplicate or non-FASTQ input')
            files.add(path)
            row[field] = path
    return expected, rows
