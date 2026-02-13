// modules/local/topology_tda.nf

process TOPOLOGY_TDA {
    tag "$meta.id"
    label 'process_medium'

    conda "environment.yml"
    container 'nf-core/manifold:dev'

    input:
    tuple val(meta), path(anndata)  // Input comes from Geometry module
    val dims                        // "0,1"

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

    topology_tda.py \\
        --input ${anndata} \\
        --output ${prefix}_topology.h5ad \\
        --homology_dimensions ${dims} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        topology_tda: \$(topology_tda.py --version | grep topology | sed 's/topology_tda.py: //g')
    END_VERSIONS
    """
}
