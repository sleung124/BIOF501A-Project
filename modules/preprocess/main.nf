params.mito_thresh = 10

process PREPROCESS_DATA {
    /* steps to preprocess visium data 
        - from visium sample, filter out all capture spots with 
          mitochondrial content above threshold
    */

    publishDir params.outputDir, mode: 'copy'

    input:
        path(script) 

    output:
        ...

    script:
        """
        Rscript $script ...
        """
}