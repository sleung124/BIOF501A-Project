/* 
    main workflow for all steps 
    TODO:
        - figure out what params are needed
*/

// import processes from specific modules
include { PREPROCESS_DATA }     from "./modules/preprocess"
// include { CELL_DECONVOLUTION }  from "./modules/cell_deconvolution"
// include { FIND_DEGS }           from "./modules/differential_expression"
// include { FIND_PATHWAYS }       from "./modules/pathway_analysis"

workflow {
    // don't know what params I want yet so not gonna put it down
    
    // path to scripts
    // preprocess_script_ch = file("$projectDir/bin/preprocess.R")
    
    // very rough outline of steps
    data = PREPROCESS_DATA(params.preprocess.MITO_THRESHOLD)
    // data_cd = CELL_DECONVOLUTION(data)
    // degs = FIND_DEGS(data_cd)
    // paths = FIND_PATHWAYS(degs)
}
