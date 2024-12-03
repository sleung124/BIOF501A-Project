process FIND_PATHWAYS {
    // steps to perform cell deconvolution on preprocessed data
    publishDir (
        path: "results/degs",
        mode: "copy",
        saveAs: { filename -> filename.equals('not_want.csv') ? null : filename }
    )
    tag "Find enriched pathways..."
    container "sleung124/spatial-pipeline:latest"
    debug "true"

    script:
        """
        pathways.r
        """
}