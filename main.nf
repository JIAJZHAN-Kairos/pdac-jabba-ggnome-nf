nextflow.enable.dsl=2

process JABBA_SAMPLE {
    tag { meta.sample }
    publishDir "${params.outdir}/samples", mode: 'copy', overwrite: false
    input:
    tuple val(meta), path(cobalt), path(sv), path(cnv), path(purity)
    path scripts
    path pins
    output:
    tuple val(meta.sample), path("result/${meta.sample}"), emit: results
    script:
    def metadata = groovy.json.JsonOutput.toJson(meta)
    """
    cat > sample.json <<'SAMPLE_METADATA'
    ${metadata}
    SAMPLE_METADATA
    run_sample.py sample.json ${scripts} ${pins} result/${meta.sample} '${cobalt}' '${sv}' '${cnv}' '${purity}'
    """
    stub:
    """
    mkdir -p result/${meta.sample}
    printf '%s\\n' '{"sample":"${meta.sample}","status":"STUB_ONLY","accepted":false}' > result/${meta.sample}/execution_status.json
    """
}

process REPORT_OUTCOMES {
    publishDir "${params.outdir}/cohort", mode: 'copy', overwrite: false
    input:
    path sample_results
    output:
    path 'cohort_outcomes.json'
    path 'cohort_outcomes.tsv'
    script:
    """
    report_outcomes.py . ${params.expected_samples}
    """
    stub:
    """
    printf '%s\\n' '{"stub":true,"accepted":0,"failed":0}' > cohort_outcomes.json
    printf 'sample\\tstatus\\nSTUB\\tSTUB_ONLY\\n' > cohort_outcomes.tsv
    """
}

process EXECUTION_STATUS_GATE {
    input:
    path summary
    output:
    path 'EXECUTION_FINISHED_REVIEW_REQUIRED'
    script:
    """
    python3 -c 'import json; d=json.load(open("${summary}")); assert d["failed"] == 0, "Scientific failures preserved in published sample results; no QC acceptance"'
    echo 'Execution finished; original full numerical and actual image QC still required' > EXECUTION_FINISHED_REVIEW_REQUIRED
    """
    stub:
    """
    echo 'STUB_ONLY; no scientific execution or acceptance' > EXECUTION_FINISHED_REVIEW_REQUIRED
    """
}

workflow {
    if (!params.input || !params.outdir || !params.container || !params.expected_samples || !params.wave_id)
        error 'Required: --input --outdir --container --expected_samples --wave_id'
    if (!(params.wave_id ==~ /[A-Za-z0-9_.-]+/)) error 'Invalid wave_id'
    rows = Channel.fromPath(params.input, checkIfExists:true).splitCsv(header:true)
        .map { row ->
            ['sample','cohort','cobalt_ratio','sv_vcf','purple_cnv','purple_purity','cobalt_sha256','sv_sha256','cnv_sha256','purity_sha256','cobalt_bytes','sv_bytes','cnv_bytes','purity_bytes'].each { k -> if (!row[k]) error "Missing ${k}" }
            if (!(row.sample ==~ /[A-Za-z0-9_.-]+/) || !(row.cohort ==~ /[A-Za-z0-9_.-]+/)) error 'Unsafe sample/cohort'
            ['cobalt','sv','cnv','purity'].each { k ->
                if (!(row["${k}_sha256"] ==~ /[a-f0-9]{64}/) || !(row["${k}_bytes"] ==~ /[0-9]+/)) error "Invalid hash/size ${row.sample}"
            }
            def meta = row.findAll { k,v -> k in ['sample','cohort'] || k.endsWith('_sha256') || k.endsWith('_bytes') }
            meta.wave_id = params.wave_id
            tuple(meta, file(row.cobalt_ratio, checkIfExists:true), file(row.sv_vcf, checkIfExists:true), file(row.purple_cnv, checkIfExists:true), file(row.purple_purity, checkIfExists:true))
        }.collect(flat:false).flatMap { xs ->
            if (xs.size() != params.expected_samples.toInteger()) error 'Unexpected samplesheet count'
            if (xs.collect { it[0].sample }.unique().size() != xs.size()) error 'Duplicate sample IDs'
            xs
        }
    JABBA_SAMPLE(rows, file("${projectDir}/assets/gadi_scripts"), file("${projectDir}/provenance/gadi_script_sha256.json"))
    result_dirs = JABBA_SAMPLE.out.results.map { sample, result -> result }.collect()
    REPORT_OUTCOMES(result_dirs)
    EXECUTION_STATUS_GATE(REPORT_OUTCOMES.out[0])
}
