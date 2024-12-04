/* 
    main workflow for all steps 
*/

// import processes from specific modules
include { PREPROCESS_DATA }     from "./modules/preprocess"
include { CELL_DECONVOLUTION }  from "./modules/cell_deconvolution"
include { FIND_DEGS }           from "./modules/differential_expression"
include { FIND_PATHWAYS }       from "./modules/pathway_analysis"

workflow {
    input_path_channel = Channel.fromPath(params.preprocess.PATH_TO_SAMPLE)
    input_filtered_feature_channel = Channel.fromPath(params.preprocess.FILTERED_FEATURE_H5)

    preprocessed = PREPROCESS_DATA(params.preprocess.MITO_THRESHOLD, params.preprocess.FILTERED_FEATURE_H5, input_path_channel)
    cell_deconv = CELL_DECONVOLUTION(params.cell_deconvolution.MAX_LDA_K, params.cell_deconvolution.RADIUS, preprocessed.processed)
    degs_found = FIND_DEGS(params.degs.QUICK_SAMPLE, params.degs.SET_SEED, params.degs.PVAL_THRESH, cell_deconv.cell_deconvolved, preprocessed.processed)
    FIND_PATHWAYS(params.pathways.SHOW_TERMS, params.pathways.NUMCHAR, params.pathways.ORDER_BY, params.pathways.SPECIES, degs_found.degs)
}
