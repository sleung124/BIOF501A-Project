// Define preprocess workflow
process PREPROCESS_DATA {
    /* steps to preprocess visium data 
        - from visium sample, filter out all capture spots with 
          mitochondrial content above threshold
    */
    publishDir (
            path: "results/preprocess_data",
            mode: "copy",
            saveAs: { filename -> filename.equals('filtered_data.rds') ? null : filename }
        )
    tag "Removing low Mitochondrial capture spots..."
    // container "sleung124/spatial-pipeline:latest"

    debug "true"
    input:
        val mito_threshold
        path path_to_sample 

    output:
        // should output the following:
        // - intermediate object for cell deconvolution (rds file)
        // - mitochondrial plot
        path "*.rds", emit: processed 
        path "*.jpg", emit: mitoplot

    script:
    // Rscript ${path} ${params.preprocess.MITO_THRESHOLD} ${params.preprocess.PATH_TO_SAMPLE}
        """
        echo this is mito threshold: "${mito_threshold}"
        preprocess.r ${mito_threshold} ${path_to_sample}
        """
}