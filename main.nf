#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 *========================================================================================
 * nf-core/manifold
 *========================================================================================
 * Pipeline for topological, geometric, and causal analysis of omics data.
 * https://github.com/nf-core/manifold
 *========================================================================================
 */

// --- HEADER LOGGING ---
log.info """\
========================================================================================
    N F - C O R E / M A N I F O L D   v 1.0 dev
========================================================================================
    Input           : ${params.input}
    Outdir          : ${params.outdir}
    Manifold Method : ${params.manifold_method}
    Causal Genes    : ${params.n_top_genes}
========================================================================================
"""

// --- MODULE IMPORTS ---
include { GEOMETRY_MANIFOLD } from './modules/local/geometry_manifold.nf'
include { TOPOLOGY_TDA      } from './modules/local/topology_tda.nf'
include { CAUSALITY_INFER   } from './modules/local/causality_infer.nf'

// --- MAIN WORKFLOW ---
workflow {

    // 1. INPUT CHANNEL INITIALIZATION
    ch_input = Channel.fromPath(params.input)
                      .map { file -> tuple([id: file.simpleName], file) }

    // 2. STEP 1: GEOMETRIC MANIFOLD LEARNING
    GEOMETRY_MANIFOLD ( 
        ch_input, 
        params.n_neighbors, 
        params.manifold_method 
    )

    // 3. STEP 2: TOPOLOGICAL DATA ANALYSIS (TDA)
    // Consumes output from Geometry module
    TOPOLOGY_TDA (
        GEOMETRY_MANIFOLD.out.h5ad,
        "0,1" 
    )

    // 4. STEP 3: CAUSAL INFERENCE
    // Consumes output from Topology module
    CAUSALITY_INFER (
        TOPOLOGY_TDA.out.h5ad,
        params.n_top_genes
    )

    // 5. OUTPUT MONITORING
    CAUSALITY_INFER.out.h5ad.view { meta, file -> 
        "Pipeline finished. Causal network inferred for: ${meta.id} -> ${file}" 
    }

}

// --- COMPLETION HANDLER ---
workflow.onComplete {
    log.info ( workflow.success ? "\nPipeline completed successfully! Results in: ${params.outdir}\n" : "\nPipeline failed.\n" )
}
