#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * -------------------------------------------------
 * nf-core/manifold: ENTRY POINT
 * -------------------------------------------------
 */

// Import the subworkflow
include { MANIFOLD_LEARNING } from './subworkflows/local/manifold_learning'

// Log info
log.info """\
         M A N I F O L D   P I P E L I N E
         ===================================
         input   : ${params.input}
         outdir  : ${params.outdir}
         methods : ${params.manifold_methods}
         """
         .stripIndent()

workflow {

    // 1. Create input channel
    // We create a "meta" map [id: 'test_sample'] to mimic nf-core standards
    ch_input = Channel.fromPath(params.input)
                      .map { file -> [ [id: 'test_sample'], file ] }

    // 2. Run the Manifold Learning subworkflow
    // We pass the input channel and the methods parameter
    MANIFOLD_LEARNING ( 
        ch_input,
        params.manifold_methods
    )

    // 3. Collect versions 
    MANIFOLD_LEARNING.out.versions.view { "Tool versions: $it" }
    
}
