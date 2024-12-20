// process for reference-free cell deconvolution
process CELL_DECONVOLUTION {
    //steps to find cell types
    publishDir (
        path: "results/cell_deconvolution",
        mode: "copy",
        saveAs: { filename -> filename.equals('deconv_results.rds') ? null : filename }
    )
    tag "Performing cell deconvolution..."
    
    input:
        val max_lda_k
        val radius
        path path_to_preprocessed

    output:
    // should output intermediate object for DEG step in the form of an RDS file
    // should output plot for cell deconvolution
        path "*.rds", emit: cell_deconvolved
        path "*.jpg", emit: deconv_plot

    script:
        """
        cell_deconvolution.r ${max_lda_k} ${radius} ${path_to_preprocessed}
        """
}