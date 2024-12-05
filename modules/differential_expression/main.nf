process FIND_DEGS{
    //steps to find DEGs from cell-abundant regions
    publishDir (
            path: "results/degs",
            mode: "copy",
            saveAs: { filename -> filename.equals('degs.rds') ? null : filename }
        )
    tag "Performing differential gene analysis..."

    input:
        val quick_sample
        val set_seed
        val p_val_thresh
        path path_to_deconv
        path path_to_preprocessed
    output:
    // should output intermediate object for pathway analysis in the form of an RDS file
    // should output volcano plot and CSV of DEGs
        path "*.rds", emit: degs
        path "*.csv"
        path "*.jpg"

    script:
        """
        find_degs.r ${quick_sample} ${set_seed} ${p_val_thresh} ${path_to_deconv} ${path_to_preprocessed}
        """
}