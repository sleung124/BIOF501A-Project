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
    
    // very rough outline of steps
    data = PREPROCESS_DATA(...)
    data_cv = CELL_DECONVOLUTION(data)
    degs = FIND_DEGS(data_cv)
    paths = FIND_PATHWAYS(degs)
}
