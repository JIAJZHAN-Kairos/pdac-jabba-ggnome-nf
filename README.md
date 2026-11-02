# PDAC JaBbA + gGnome Nextflow

A portable wrapper around the pinned, previously executed Gadi scientific scripts. The model, purity/ploidy inputs, CBS adaptation, native event classifier and scientific source files are preserved byte-for-byte. This is a custom execution wrapper because the existing implementation is a PBS/Singularity workflow and the live workflow catalog has no equivalent JaBbA/gGnome entry.

## Inputs

One CSV with columns `sample,cohort,cobalt_ratio,sv_vcf,purple_cnv,purple_purity,cobalt_sha256,sv_sha256,cnv_sha256,purity_sha256,cobalt_bytes,sv_bytes,cnv_bytes,purity_bytes`. Four original COBALT/ESVEE/PURPLE files per tumor; no BAM/FASTQ. GRCh38 coordinates and original SHA256 values are mandatory. Purity, ploidy and gender are read from the original single-row PURPLE purity file.

## Execution

```bash
nextflow run JIAJZHAN-Kairos/pdac-jabba-ggnome-nf -r <immutable-commit> -profile aws \
  --input s3://<private-bucket>/<prefix>/samplesheet.csv \
  --outdir s3://<private-bucket>/<unique-run>/results \
  --container <private-ecr-image>@sha256:<digest> \
  --expected_samples 2 --wave_id migration-pilot-20261030
```

Use the user-selected Seqera AWS Batch compute environment. Without Fusion, AWS CLI stages S3 inputs and outputs. The job asks for four CPUs, 16 GiB, eight hours, four solver threads and CPLEX max.mem 12 GB; tilim 6000, epgap 1e-4, slack penalty 100 are fixed in the pinned R entrypoint. Maximum 20 concurrent sample processes. No automatic resource escalation or scientific retry.

## Runtime

The private image uses the exact Docker base digest from the original Gadi SIF and the unchanged installed R libraries/CPLEX runtime. The runtime bundles and proprietary binaries must stay outside Git and public container registries. `containers/Dockerfile` expects the hash-verified private `runtime_pinned.tar.gz` plus `cplex_include_lib_original.tar.gz`, which together include the original CPLEX `bin`, `include`, and `lib` directories. Their SHA256 values are recorded in `provenance/runtime_build.json`. CPLEX use remains subject to the existing entitlement. Record the immutable ECR image digest for every run.

The first migration image contained the CPLEX executable but omitted its headers and static libraries. JaBbA rejected that installation before optimization. Version 0.1.1 checks the required directories and original core-file SHA256 values before processing inputs; empty directories do not satisfy this check. Existing failed run outputs remain untouched. Stage logs are retained in the result and streamed to private AWS logs for diagnosis.

## Failure preservation and QC

Each sample records its original worker status, logs, original graph/checkpoint files, scientific-script hashes and output hashes. A sample process stages out these files even if its scientific subprocess fails; the final execution-status gate then fails the pipeline when any sample failed. This prevents lost original diagnostics without presenting failed inference as a successful analysis.

Numerical gate: every final component status 101/102, finite gap <=1e-4, nonnegative integer known CN, explicit unknown bp, native internal balance exactly zero. It does not waive strict numeric checks because upstream export uses a tolerance. Actual review of CN, highest CN, every native complex-event PNG, original coverage support and unresolved boundaries is still required. The pipeline never writes ACCEPTED or a pilot-QC pass marker. `OUTPUTS_COMPLETE_REVIEW_PENDING` is not acceptance.

PURPLE comparison is consistency, not independent truth. Native dm/cpxdm is not validated ecDNA; BFB/groupings are candidate mechanisms and native chromothripsis is exploratory rather than ShatterSeek. Unknown/failed/uncalled is not a biological negative. Preserve all original architecture and regional-CN limitations.

## Sample ownership

A samplesheet is not a reservation. Before any cloud submission, reconcile actual Gadi receipts and install an exclusive migration ownership record under the same controller lock. Never rerun accepted results or cloud-submit a sample already assigned to a healthy Gadi job. Existing recovery attempts remain separate and retain their raw/checkpoint provenance.

## Provenance

`provenance/gadi_script_sha256.json` is checked before execution. The Gadi originals and its accepted graphs are never modified by this repository. Stub runs test orchestration only and write STUB_ONLY with accepted=false.

## Explicit recovery wave

Version 0.2.0 adds `--recovery true` for a new, separately owned samplesheet. It adds `recovery_profile`, `recovery_spec` and `recovery_spec_sha256` columns to the original input contract. A private per-sample JSON spec binds all four original input hashes, mode, resource profile, adapter hash and any original checkpoint hashes. Existing completed samples must not enter this manifest.

The `standard` profile remains 4 CPU/16 GiB/8 h with CPLEX 12 GB/tilim 6000. The `recovery` profile uses 12 CPU/64 GiB/48 h with CPLEX 48 GB/tilim 21600, and `hard` uses 24 CPU/96 GiB/48 h with CPLEX 72 GB/tilim 43200. Normal recovery and hard recovery have concurrency caps of three and two; the global executor queue size remains 20. Select an existing On-Demand queue for long recovery processes in the launch configuration; do not change a shared compute environment or its other jobs.

`NATIVE_WARMSTART` requires the complete original `jabba.raw.rds`, `karyograph.rds`, `edge.nudge.rds`, `ab.force.rds` and `ab.exclude.rds`, verified against their original input identity and SHA256. It uses the same native `ramip_stub` then `JaBbA(overwrite=FALSE)` method as the existing Gadi recovery. A preprocessing mipstart or component solution is never a substitute. `LONG_FRESH` starts an independent native fit when no complete checkpoint exists. No graph splicing, changed purity/ploidy, relaxed gap/status, or altered event classifier is used.

Recovery infrastructure failures do not cancel unrelated samples. The cohort report explicitly lists any sample with no published worker record as requiring scheduler review. The final gate still fails when any requested sample lacks complete output; successful execution remains separate from manual numerical/image/coverage acceptance.

The portable Python QC environment uses Python 3.12.8 rather than the Gadi controller Python 3.13.9. NumPy 1.26.4, pandas 2.2.2 and matplotlib 3.9.2 match the observed Gadi versions. R model preparation and inference use the unchanged original R4.0.5 stack. This portability difference is recorded in provenance/runtime_build.json and requires actual migration-pilot validation.
