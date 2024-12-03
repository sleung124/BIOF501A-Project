// Define preprocess workflow
process PREPROCESS_DATA {
    /* steps to preprocess visium data 
        - from visium sample, filter out all capture spots with 
          mitochondrial content above threshold
    */
    publishDir (
            path: "results/preprocess_data",
            mode: "copy",
            saveAs: { filename -> filename.equals('not_want.csv') ? null : filename }
        )
    tag "Removing low Mitochondrial capture spots..."
    container "sleung124/spatial-pipeline:latest"
    // publishDir params.outputDir, mode: 'copy'
    // path "${workflow.projectDir}/bin/preprocess.r"
    debug "true"
    input:
        // path(script)
        val mito_threshold
        path path_to_sample 

    output:
        // path("dummy.csv")
        // path("not_want.csv")
        // stdout
        // should output the following:
        // - intermediate object for cell deconvolution (rds file)
        // - mitoplot
        // path "*.rds", emit: processed 
        // path "*.jpg", emit: mitoplot
        stdout

    script:
    // Rscript ${path} ${params.preprocess.MITO_THRESHOLD} ${params.preprocess.PATH_TO_SAMPLE}
        """
        echo this is mito threshold: "${mito_threshold}"
        preprocess.r ${mito_threshold} ${path_to_sample}
        """
}