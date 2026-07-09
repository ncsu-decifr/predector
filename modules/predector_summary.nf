process PREDICTOR_SUMMARY {

    tag "$name"

    publishDir "${params.outdir}/${name}/summary", mode: 'copy'

    input:
    tuple val(name), path(tsv_files)

    output:
    tuple val(name), path("summary/*")

    script:
    """
    mkdir -p summary

    Rscript ${projectDir}/bin/predector_summary.R \
        --indir . \
        --outdir summary \
        --assembly "${params.summary_assembly}"
    """
}