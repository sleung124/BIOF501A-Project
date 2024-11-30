// Define preprocess workflow

process PREPROCESS_DATA {
    /* steps to preprocess visium data 
        - from visium sample, filter out all capture spots with 
          mitochondrial content above threshold
    */
    tag "Removing low Mitochondrial capture spots..."
    container "sleung124/spatial-pipeline:latest"
    // publishDir params.outputDir, mode: 'copy'
    // path "${workflow.projectDir}/bin/preprocess.r"
    debug "true"
    input:
        // path(script)
        val(mito_threshold) 

    output:
        val("hello")
        // should output the following:
        // - intermediate object for cell deconvolution (rds file)
        // - mitoplot
        // path("*.rds"), emit processed 
        // path("${publishDir}/preprocess/*.jpg")

    script:
    // Rscript ${path} ${params.preprocess.MITO_THRESHOLD} ${params.preprocess.PATH_TO_SAMPLE}
        """
        echo this is mito threshold: "${mito_threshold}"
        temp.r
        """
}