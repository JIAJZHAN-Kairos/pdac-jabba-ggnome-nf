process RECOVERY_SAMPLE {
    tag { meta.sample }
    publishDir "${params.outdir}/samples", mode: 'copy', overwrite: false
    errorStrategy 'ignore'
    input:
    tuple val(meta), path(cobalt), path(sv), path(cnv), path(purity), path(spec)
    path scripts
    path pins
    path cplex_pins
    path adapter
    output:
    tuple val(meta.sample), path("result/${meta.sample}"), emit: results
    script:
    def metadata = groovy.json.JsonOutput.toJson(meta)
    """
    export OMP_NUM_THREADS=${task.cpus}
    cat > sample.json <<'SAMPLE_METADATA'
    ${metadata}
    SAMPLE_METADATA
    run_sample.py sample.json ${scripts} ${pins} ${cplex_pins} result/${meta.sample} '${cobalt}' '${sv}' '${cnv}' '${purity}' '${spec}' '${adapter}'
    """
    stub:
    """
    mkdir -p result/${meta.sample}
    printf '%s\\n' '{"sample":"${meta.sample}","status":"STUB_ONLY","accepted":false,"profile":"${meta.recovery_profile}"}' > result/${meta.sample}/execution_status.json
    """
}
