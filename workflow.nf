/* 
    main workflow for all steps 
    TODO:
        - figure out what params are needed
*/


// preprocess step 
process PREPROCESS_DATA {
    // steps to preprocess visium data 
}

// cell deconvolution step
process CELL_DECONVOLUTION {
    // steps to perform cell deconvolution on preprocessed data
}

// differential gene analysis step
process FIND_DEGS{
    // steps to find DEGs between cell enriched regions
}

// patway analysis step
process FIND_PATHWAYS {
    //steps to find enriched pathways from DEGs
}

workflow {
    // don't know what params I want yet so not gonna put it down
    
    // very rough outline of steps
    data = PREPROCESS_DATA(...)
    data.cv = CELL_DECONVOLUTION(data)
    degs = FIND_DEGS(data.cv)
    paths = FIND_PATHWAYS(degs)
}
