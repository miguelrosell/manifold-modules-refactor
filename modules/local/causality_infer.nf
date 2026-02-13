// modules/local/causality_infer.nf

process CAUSALITY_INFER {
    tag "$meta.id"
    label 'process_medium'

    conda "environment.yml"
    container 'nf-core/manifold:dev'

    input:
    tuple val(meta), path(anndata)  // Input comes from Topology module
    val n_top_genes                 // Parameter for network size

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    export NUMBA_CACHE_DIR=/tmp

    causality_infer.py \\
        --input ${anndata} \\
        --output ${prefix}_causal.h5ad \\
        --n_top_genes ${n_top_genes} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        causality_infer: \$(causality_infer.py --version | grep causality | sed 's/causality_infer.py: //g')
    END_VERSIONS
    """
}
