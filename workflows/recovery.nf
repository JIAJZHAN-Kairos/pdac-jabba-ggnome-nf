include { RECOVERY_SAMPLE as FIRST_RUN; RECOVERY_SAMPLE as NATIVE_RECOVERY; RECOVERY_SAMPLE as HARD_RECOVERY } from './recovery_sample'

process RECOVERY_REPORT {
    publishDir "${params.outdir}/cohort", mode: 'copy', overwrite: false
    input:
    path sample_results
    path manifest
    output:
    path 'cohort_outcomes.json'
    path 'cohort_outcomes.tsv'
    script:
    """
    report_recovery_outcomes.py . '${manifest}' ${params.expected_samples}
    """
    stub:
    """
    printf '%s\\n' '{"stub":true,"accepted":0,"failed":0}' > cohort_outcomes.json
    printf 'sample\\tstatus\\nSTUB\\tSTUB_ONLY\\n' > cohort_outcomes.tsv
    """
}

process RECOVERY_GATE {
    input:
    path summary
    output:
    path 'EXECUTION_FINISHED_REVIEW_REQUIRED'
    script:
    """
    python3 -c 'import json; d=json.load(open("${summary}")); assert d["failed"] == 0, "Failures preserved in cohort report; no QC acceptance"'
    echo 'Execution finished; original numerical/image/coverage review still required' > EXECUTION_FINISHED_REVIEW_REQUIRED
    """
    stub:
    """
    echo 'STUB_ONLY; no scientific execution or acceptance' > EXECUTION_FINISHED_REVIEW_REQUIRED
    """
}

workflow RECOVERY76 {
    rows = Channel.fromPath(params.input, checkIfExists:true).splitCsv(header:true)
        .map { row ->
            ['sample','cohort','cobalt_ratio','sv_vcf','purple_cnv','purple_purity','cobalt_sha256','sv_sha256','cnv_sha256','purity_sha256','cobalt_bytes','sv_bytes','cnv_bytes','purity_bytes','recovery_spec','recovery_spec_sha256','recovery_profile'].each { k -> if (!row[k]) error "Missing ${k}" }
            if (!(row.sample ==~ /[A-Za-z0-9_.-]+/) || !(row.cohort ==~ /[A-Za-z0-9_.-]+/)) error 'Unsafe sample/cohort'
            if (!(row.recovery_profile in ['standard','recovery','hard'])) error 'Invalid recovery profile'
            ['cobalt','sv','cnv','purity'].each { k ->
                if (!(row["${k}_sha256"] ==~ /[a-f0-9]{64}/) || !(row["${k}_bytes"] ==~ /[0-9]+/)) error "Invalid hash/size ${row.sample}"
            }
            if (!(row.recovery_spec_sha256 ==~ /[a-f0-9]{64}/)) error 'Invalid recovery spec hash'
            def meta = row.findAll { k,v -> k in ['sample','cohort','recovery_profile','recovery_spec_sha256'] || k.endsWith('_sha256') || k.endsWith('_bytes') }
            meta.wave_id = params.wave_id
            tuple(meta, file(row.cobalt_ratio, checkIfExists:true), file(row.sv_vcf, checkIfExists:true), file(row.purple_cnv, checkIfExists:true), file(row.purple_purity, checkIfExists:true), file(row.recovery_spec, checkIfExists:true))
        }.collect(flat:false).flatMap { xs ->
            if (xs.size() != params.expected_samples.toInteger()) error 'Unexpected samplesheet count'
            if (xs.collect { it[0].sample }.unique().size() != xs.size()) error 'Duplicate sample IDs'
            xs
        }.branch { first: it[0].recovery_profile == 'standard'; normal: it[0].recovery_profile == 'recovery'; hard: it[0].recovery_profile == 'hard' }
    scripts = file("${projectDir}/assets/gadi_scripts")
    pins = file("${projectDir}/provenance/gadi_script_sha256.json")
    cplex = file("${projectDir}/provenance/cplex_core_sha256.json")
    adapter = file("${projectDir}/assets/recovery/run_sample_recovery.R")
    FIRST_RUN(rows.first, scripts, pins, cplex, adapter)
    NATIVE_RECOVERY(rows.normal, scripts, pins, cplex, adapter)
    HARD_RECOVERY(rows.hard, scripts, pins, cplex, adapter)
    results = FIRST_RUN.out.results.mix(NATIVE_RECOVERY.out.results, HARD_RECOVERY.out.results).map { sample, result -> result }.collect()
    RECOVERY_REPORT(results, file(params.input))
    RECOVERY_GATE(RECOVERY_REPORT.out[0])
}
