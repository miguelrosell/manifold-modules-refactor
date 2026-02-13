// modules/local/geometry_manifold.nf

process GEOMETRY_MANIFOLD {
    tag "$meta.id"
    label 'process_medium'

    // Define the container to use. 
    // In the future, this will point to quay.io/nf-core/manifold, but for now we use local context.
    conda "environment.yml"
    container 'nf-core/manifold:dev'

    input:
    tuple val(meta), path(anndata)  // The input file (e.g., sample1.h5ad)
    val n_neighbors                 // Parameter for graph construction
    val method                      // 'phate' or 'diffmap'

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad     // The processed file
    path "versions.yml"            , emit: versions // Version tracking file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    export NUMBA_CACHE_DIR=/tmp

    # Run the Python script located in bin/
    geometry_manifold.py \\
        --input ${anndata} \\
        --output ${prefix}_geometric.h5ad \\
        --n_neighbors ${n_neighbors} \\
        --method ${method} \\
        $args

    # Capture versions for provenance
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        geometry_manifold: \$(geometry_manifold.py --version | grep geometry | sed 's/geometry_manifold.py: //g')
    END_VERSIONS
    """
}
