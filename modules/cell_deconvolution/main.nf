process CELL_DECONVOLUTION {
    //steps to find cell types
    publishDir (
        path: "results/cell_deconvolution",
        mode: "copy",
        saveAs: { filename -> filename.equals('not_want.csv') ? null : filename }
    )
    tag "Performing cell deconvolution..."
    container "sleung124/spatial-pipeline:latest"
    debug "true"

    script:
        """
        cell_deconvolution.r
        """
}