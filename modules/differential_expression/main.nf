process FIND_DEGS{
    //steps to find enriched pathways from DEGs
    publishDir (
            path: "results/degs",
            mode: "copy",
            saveAs: { filename -> filename.equals('not_want.csv') ? null : filename }
        )
    tag "Performing differential gene analysis..."
    container "sleung124/spatial-pipeline:latest"
    debug "true"

    script:
        """
        find_degs.r
        """
}