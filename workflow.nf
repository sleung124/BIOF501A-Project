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
    PREPROCESS_DATA(params.preprocess.MITO_THRESHOLD, input_path_channel)
    CELL_DECONVOLUTION()
    FIND_DEGS()
    FIND_PATHWAYS()
}
