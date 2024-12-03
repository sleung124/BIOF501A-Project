process FIND_PATHWAYS {
    // steps to perform cell deconvolution on preprocessed data
    publishDir (
        path: "results/pathways",
        mode: "copy"
    )
    tag "Find enriched pathways..."
    // container "sleung124/spatial-pipeline:latest"
    debug "true"
    input:
        val pval_thresh 
        val show_terms 
        val numchar 
        val order_by 
        val dbs 
        path path_to_degs
    output:
        path "*.jpg"

    script:
        """
        pathways.r ${pval_thresh} ${show_terms} ${numchar} ${order_by} ${dbs} ${path_to_degs}
        """
}