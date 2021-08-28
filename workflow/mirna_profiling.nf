nextflow.enable.dsl=2


process alignment_stat {

  publishDir "$params.out_dir", mode: 'copy', overwrite: true

  input:
  tuple val(sampleId), path(bam_plus)

  output:
  tuple val(sampleId), path ("$sampleId") 

  script:

  def config =  params.config_file != null ? "--config_file ${projectDir}/$params.config_file" : ""
  def sampleBam = bam_plus.find{it != null && it.name.endsWith(".bam") }

  """
  alignment_stat --bam $sampleBam --ucsc $params.ucsc --mirbase $params.mirbase --species $params.species --nd --out_dir $sampleId $config
  """
}


process expression_matrix {

  publishDir "$params.out_dir", mode: 'copy', overwrite: true

  input: 
  tuple val(sampleId), path (sample_dir) 
  
  output:
  path (sample_dir) 

  script:

  def config =  params.config_file != null ? "--config_file ${projectDir}/$params.config_file" : ""
  
  """
  expression_matrix --proj_dir $sample_dir --mirbase $params.mirbase --species $params.species --out_dir $sample_dir $config
  """
}

process project_summary {

  publishDir "$params.out_dir", mode: 'copy', overwrite: true

  input:
  val results 

  script:

  def proj_name = params.proj_name == null ? new File("$params.out_dir").name : params.proj_name
 
  """
  project_summary --proj_dir $params.runtime_out_dir --name $proj_name
  """
}

