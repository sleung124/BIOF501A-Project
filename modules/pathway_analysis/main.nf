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
        val show_terms 
        val numchar 
        val order_by 
        val species
        path path_to_degs
    output:
        path "*.jpg"

    script:
        """
        pathways.r  ${show_terms} ${numchar} ${order_by} ${species} ${path_to_degs}
        """
}