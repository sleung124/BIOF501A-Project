process FIND_DEGS{
    //steps to find enriched pathways from DEGs
    publishDir (
            path: "results/degs",
            mode: "copy",
            saveAs: { filename -> filename.equals('degs.rds') ? null : filename }
        )
    tag "Performing differential gene analysis..."
    // container "sleung124/spatial-pipeline:latest"
    debug "true"
    input:
        val quick_sample
        path path_to_deconv
        path path_to_preprocessed
    output:
        path "*.rds", emit: degs
        path "*.csv"

    script:
        """
        find_degs.r ${quick_sample} ${path_to_deconv} ${path_to_preprocessed}
        """
}