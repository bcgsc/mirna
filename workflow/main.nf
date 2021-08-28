/*
  TODO: Replace container path with published container URI
*/

nextflow.enable.dsl=2

include { alignment_stat ; expression_matrix ; project_summary } from './mirna_profiling.nf'

def helpMessage() {
  log.info """
        Usage:
        nextflow run main.nf -profile cluster,container --proj_dir /path/bam_dir --out_dir /path/project_dir
        
        nextflow run main.nf -profile local -params-file params.yaml

        Mandatory arguments:
         -profile                   cluster,container : run container on cluster
                                    container :  run container on localhost
                                    cluster,local : run installed package (venv) on cluster
                                    local :  run installed package (venv) on localhost

         --proj_dir                 Directory where bams are stored 
         --out_dir                  Output directory for profiling results
         --ucsc                     UCSC genome or database name. Default: hg38
         --mirbase                  miRBase release or database name. Default: mirna_21
         --species                  miRBase species code.  Default: hsa

       Optional arguments:
        --sample_csv                CSCV file with sample, bam and adapter columns
        --proj_name                 Project name to use in graphs
        --config_file               Custom config yaml for profiling. Must be in the same directory as main.nf, 
                                    if run with container profile

        """
}

if (params.help) {
    helpMessage()
    exit 0
}

workflow {
    if ( params.sample_csv != null ) {
        data = Channel
        .fromPath(params.sample_csv)
        .splitCsv(header:true)
        .map{ row-> row.adapter == null ? tuple(row.sample, [file(row.bam)]) : tuple(row.sample, [file(row.bam), file(row.adapter)]) }
    } else {
        data = Channel
        .fromPath(["$params.proj_dir/**.bam", "$params.proj_dir/**adapter*report*"])
        .map { file ->
        def key = file.name.replaceFirst(/[\._].+$/, '')
        return tuple(key, file)
        }
        .groupTuple()
    }

    expression_matrix(alignment_stat(data))
    project_summary(expression_matrix.out.toList())
}
   

