/* 
    main workflow for all steps 
    TODO:
        - figure out what params are needed
*/

// import processes from specific modules
include { PREPROCESS_DATA }     from "./modules/preprocess"
include { CELL_DECONVOLUTION }  from "./modules/cell_deconvolution"
include { FIND_DEGS }           from "./modules/differential_expression"
include { FIND_PATHWAYS }       from "./modules/pathway_analysis"

workflow {
    // don't know what params I want yet so not gonna put it down
    input_path_channel = Channel.fromPath(params.preprocess.PATH_TO_SAMPLE)
    // very rough outline of steps
    preprocessed = PREPROCESS_DATA(params.preprocess.MITO_THRESHOLD, input_path_channel)
    cell_deconv = CELL_DECONVOLUTION(params.cell_deconvolution.MAX_LDA_K, params.cell_deconvolution.RADIUS, preprocessed.processed)
    degs = FIND_DEGS(params.degs.QUICK_SAMPLE, cell_deconv.cell_deconvolved, preprocessed.processed)
    FIND_PATHWAYS(params.pathways.PVAL_THRESH, params.pathways.SHOW_TERMS, params.pathways.NUMCHAR, params.pathways.ORDER_BY, params.pathways.DBS, degs.degs)
}
