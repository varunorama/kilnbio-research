#!/usr/bin/env nextflow
// Native KilnBio study-level dispatch. Scientific work occurs in direct modules.
nextflow.enable.dsl=2
params.manifest=null
params.outdir=null
params.max_parallel=1

process ANALYZE {
    tag "$id:$mode"
    maxForks params.max_parallel as int
    stageInMode 'symlink'
    publishDir "${params.outdir}", mode:'copy', overwrite:false
    input:
    tuple val(id),val(mode),path(config)
    output:
    tuple val(id),path("${id}")
    script:
    def root=projectDir.parent.toString().replace("'", "'\\''")
    """
    bash '${root}/run.sh' '${mode}' '${config}' '${id}'
    """
}

workflow {
    if(!params.manifest || !params.outdir)error 'Use --manifest jobs.csv --outdir NEW_OUTPUT'
    if(file(params.outdir).exists())error 'Output exists; choose a new output directory'
    def manifest=file(params.manifest,checkIfExists:true)
    def seen=[]
    jobs=Channel.fromPath(manifest).splitCsv(header:true).map { row ->
        if(!(row.id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || seen.contains(row.id))error 'Unique safe job IDs required'
        seen.add(row.id)
        if(!(row.mode in ['de','dm','featurecounts','enrich','scrna','pseudobulk','import-mods','ont']))error 'Unsupported analysis mode'
        def config=file(manifest.parent.resolve(row.config),checkIfExists:true)
        if(config.name.contains("'"))error 'Apostrophes in config names are unsupported'
        tuple(row.id,row.mode,config)
    }
    ANALYZE(jobs)
}
