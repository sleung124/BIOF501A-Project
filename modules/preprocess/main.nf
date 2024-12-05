// process for filtering out low-quality capture spots
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

    input:
        val mito_threshold
        val filtered_feature_h5 
        path path_to_sample 

    output:
        // should output the following:
        // - intermediate object for cell deconvolution (rds file)
        // - mitochondrial plot
        path "*.rds", emit: processed 
        path "*.jpg", emit: mitoplot

    script:
        """
        echo this is mito threshold: "${mito_threshold}"
        preprocess.r ${mito_threshold} ${filtered_feature_h5} ${path_to_sample}
        """
}