# SpliceAI search-distance experiment

This experiment compares `D=50`, `4999`, `5000`, and `10000` using the exact
SpliceAI container and reference-transcript annotation pinned by ARIANE. It is
an isolated research check. It does not change the Appendix J production profile.

## Method

`scripts/compare_spliceai_distances.py` calls the installed upstream SpliceAI
scorer directly, without an HTTP server or result cache. It captures each of
the five ensemble models' REF and ALT outputs, averages them in the same order
as the scorer, and maps negative-strand outputs back to genomic coordinates.
Every reconstructed DS, DP, REF and ALT extremum must match the upstream
formatted output. The two existing service validation variants must also match
their versioned D=10000 fixtures.

The cases include five repository-backed SNVs and one deliberately synthetic
SNV at the midpoint of BRCA2's longest intron in the pinned annotation. The
synthetic variant's position and reference base are obtained from the actual
annotation and FASTA. This is not a claim that the synthetic variant occurs in
patients. Every reference allele is checked before inference.

For each case and distance, the experiment saves full unrounded REF/ALT arrays
and records:

- all four delta scores and their argmax positions;
- REF and ALT probabilities at each maximum;
- differences between overlapping prediction windows;
- the exact model-boundary positions at +/-4999, +/-5000 and +/-5001;
- maximum absolute delta and exact/tolerance-based nonzero counts beyond 5000;
- raw and three-decimal SpliceAI bands at 0.1 and 0.2, without assigning criteria;
- elapsed time after one untimed warmup per input shape.

The default run has two repeats, with cyclic and reversed distance ordering.
Numerical comparisons use absolute tolerance 1e-6, and separately record
three-decimal equality. A changed DP can result from an argmax tie or tiny
prediction differences; it is not by itself evidence of a changed splice event.
Timing describes this machine and configuration, not production API latency.

The report includes image digest, model commit, model-file checksums, FASTA and
annotation checksums, script checksum, package versions, threads, and checksums
of the retained raw arrays. The output directory must be new. Failed runs and
baseline mismatches are explicit, rather than accepted as equivalent results.

## Results: 2026-09-15

The [machine-readable report](../../../data/spliceai/distance_comparison_20260915.json)
records a completed local run: six SNVs, four distances, two repeats (48 fresh
timed evaluations), plus four excluded warmups. Both existing D=10000 validation
fixtures passed for all their delta, REF and ALT component scores.

Maximum unrounded delta scores, displayed here to nine decimal places:

| Variant | D=50 | D=4999 | D=5000 | D=10000 |
| --- | ---: | ---: | ---: | ---: |
| BRCA1 c.4185G>A | 0.929783702 | 0.929783702 | 0.929783702 | 0.929783702 |
| BRCA2 c.7805+9T>G | 0.049350619 | 0.206400454 | 0.206400454 | 0.206400454 |
| BRCA1 c.134+10A>C | 0.000443757 | 0.002829075 | 0.002829075 | 0.002829075 |
| BRCA1 c.135-13T>G | 0.202782929 | 0.233124763 | 0.233124763 | 0.233124763 |
| BRCA2 c.1909+1G>A | 0.991917431 | 0.991917431 | 0.991917431 | 0.991917431 |
| Synthetic BRCA2 chr13:32387417 G>A | 0.000162296 | 0.000295192 | 0.000295192 | 0.000295192 |

The synthetic position lies midway through MANE intron 24, interval
chr13:32380146-32394688 (1-based inclusive), length 14,543 bp.

For comparisons among D=4999, D=5000 and D=10000:

- All six maximum delta scores were exactly equal at the captured precision.
- All four component scores agreed after service rounding, and all four DP
  positions were unchanged. The largest unrounded component-score difference
  was 8.731149137e-11, in the synthetic variant's donor-loss component.
- Overlapping raw REF/ALT probabilities were not universally bit-identical:
  the largest absolute difference was 1.192092896e-7. This is below the chosen
  tolerance and did not change any maximum score or threshold band.
- Every donor/acceptor delta beyond +/-5000 was exactly zero in the captured
  ensemble arrays. The two exact boundary positions at +/-5000 also had zero
  delta for these cases. This is an observation for this pilot, not a universal
  assertion for all possible variants or numerical implementations.
- Repeated evaluations at the same distance produced identical arrays.

The smaller-window control demonstrably changed results. For
BRCA2 c.7805+9T>G, donor gain at offset -109 has delta 0.206400454. D=50 excludes
that site and yields maximum 0.049350619 from donor loss at -9 instead. Thus the
test detects a real search-window effect, including a change of SpliceAI band.
No ARIANE evidence criteria or classifications were changed by this experiment.

Median warm elapsed time over 12 evaluations per distance was 1.21 s (D=50),
1.59 s (D=4999), 1.65 s (D=5000), and 2.20 s (D=10000). In this run the ratio
D10000/D4999 was 1.380, approximately 38% slower. These are descriptive local
CPU timings with a four-CPU container limit and four TensorFlow intra-op threads;
they are not a measured production-server slowdown or a general speed guarantee.

The original report and 48 compressed raw-array files remain locally under
`data/spliceai/build/distance-20260915/` (about 8.3 MiB in total). The versioned
report is an exact copy of that run's `report.json`; its array filenames refer
to the original experiment directory. Recompute using the procedure below to
obtain another complete raw dataset.

Conclusion: this pilot found no additional maximum delta signal from expanding
D=4999 to D=10000, including at a tested position in a long intron. It supports
practical agreement for the tested SNVs, with measurable additional computation,
but does not prove universal equality or justify changing the binding Appendix J
profile without resolving the methodological question.

## Reproduce

From the repository root on Linux, with Docker running and the pinned image
available:

```bash
mkdir -p data/spliceai/build
image='docker.io/weisburd/spliceai-38@sha256:c2bb3d5c65eee01087b5dd130d72233fc812d158dce9779ae63117e826ec8399'
docker run --rm --init --network none \
  --cap-drop ALL --security-opt no-new-privileges:true \
  --cpus 4 --memory 6g \
  --env PYTHONDONTWRITEBYTECODE=1 --env OMP_NUM_THREADS=4 \
  --env TF_NUM_INTRAOP_THREADS=4 --env TF_NUM_INTEROP_THREADS=1 \
  --mount "type=bind,src=$PWD,dst=/workspace,readonly" \
  --mount "type=bind,src=$PWD/data/spliceai/build,dst=/experiment-output" \
  --entrypoint python "$image" \
  /workspace/scripts/compare_spliceai_distances.py \
  --fasta /hg38.fa.gz \
  --output /experiment-output/distance-new-run \
  --image-digest "$image" --repeats 2 --threads 4
```

Choose a fresh output directory for each run. The container has no network and
the repository mount is read-only. Only the separate experiment output mount
is writable. No running ARIANE service is contacted.

Analysis-helper tests use artificial arrays, including a deliberately introduced
signal at exactly +5000 that D=4999 must miss. They verify the analysis itself;
they are separate from the real-model experiment:

```bash
python -m pytest tests/test_spliceai_distance_experiment.py -q
```

## Interpretation limits

This pilot can demonstrate an observed difference or agreement for its six SNVs.
It cannot prove universal equivalence for all SNVs, cover indels, establish a
clinical outcome, or determine what the ENIGMA authors intended in Appendix J.
In particular, model context and scored search radius remain distinct settings.

## Transcript-wide, resumable validation

The six-SNV run above is a pilot, not a representative statistical sample.
The new `--cohort stratified` mode builds a reproducible engineering-validation
cohort from the pinned FASTA and full reference transcript annotations:

- Three uniformly sampled canonical positions in each of all 98 exon/intron
  regions: BRCA1 23 exons + 22 introns, BRCA2 27 exons + 26 introns. Every
  sampled position gets all three alternative SNV alleles. Exonic UTR sequence
  is included; this is a model experiment, not ARIANE classification coverage.
- Targeted splice-boundary positions and offsets around 4999/5000/5001 bp.
- Simple insertions and deletions of 1, 2, 3 and 6 bases on both genes/strands.
  Complex delins/MNVs are explicitly excluded; indel REF/ALT arrays are aligned
  using the pinned upstream algorithm before comparison.
- The two pinned validation controls, including D=50 as a positive control.

With the default seed and three random positions per region, the generated plan
contains **4,047 distinct alleles at 1,315 positions**: 3,919 SNVs, 64 deletions
and 64 insertions. Of these, 882 SNVs cover 294 randomly drawn positions; the
rest are targeted stress cases or controls, with overlap/deduplication recorded.
The manifest reports each region's coverage and any unavailable reference bases.

All cases compare D=4999, 5000 and 10000. Only controls also use D=50, and only
the first 12 cases repeat for timing/repeatability. This requires 12,181 scored
evaluations plus warm-up, roughly seven hours at two seconds per evaluation.
This is a planning estimate, not a measured full-run duration. Increasing
`--positions-per-region` expands random coverage without changing the design.

### Run on Windows without supervising the conversation

Docker must be running, with the pinned image available. From the repository:

```powershell
# Prepare and inspect plan.json, without loading the models.
.\scripts\run_spliceai_distance_scan.ps1 -Action Plan

# Run independently in the background, bounded to four CPUs and 6 GiB RAM.
.\scripts\run_spliceai_distance_scan.ps1 -Action Run -Detached

# Read the compact progress summary, without contacting the model.
.\scripts\run_spliceai_distance_scan.ps1 -Action Status

# After interruption, continue only unfinished cases.
.\scripts\run_spliceai_distance_scan.ps1 -Action Resume -Detached
```

If Windows blocks unsigned local scripts, invoke this file in a process-scoped
PowerShell session, without changing the machine's execution policy:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\run_spliceai_distance_scan.ps1 -Action Run -Detached
```

For a bounded test, use `-Action Run -MaxCases 6` or `-MaxSeconds 300`, then
`-Action Resume`. Budgets are checked between complete cases; model startup
and warm-up are additional. To request more random loci, use e.g.
`-PositionsPerRegion 10 -OutputName distance-denser` consistently for Plan,
Run and Resume. The Linux Docker command above supports the same mode by adding
`--cohort stratified`, with optional `--prepare-only`, `--resume`, `--max-cases`
and `--max-seconds`. Use a new directory for a changed design.

The calculation uses no LLM calls, public API, production ARIANE instance or
classification cache. Only a progress line per 100 completed cases is printed.
Results are stored under `data/spliceai/build/<OutputName>/`:

- `plan.json`: full cohort, seed, settings, reference/model/source checksums.
- `results.sqlite3`: per-case scores, differences, timing, failures and provenance;
  committed after each complete case. Interrupted work resumes without repeating
  successfully committed cases. Source/model/configuration changes reject resume.
- `summary.json` and `summary.md`: successful, failed and unfinished counts,
  differences, timing and coverage by gene/type/region. These are periodic
  snapshots; SQLite remains authoritative if the process is killed abruptly.
- `arrays/`: raw model outputs for changed/anomalous cases by default, not
  thousands of unnecessary files. `--save-arrays all` retains every case.

Errors never count as agreement or as a zero score. Retrying failed cases is
explicit (`-Action Resume -RetryErrors`). A control-score mismatch stops the run.
Successful records cannot be overwritten with a different result. An exclusive
lock prevents two runners from sharing the same experiment directory.

Interpret differences separately for random loci, targeted cases, SNVs and
indels. Three ALT alleles at one locus are correlated; repeats are not additional
samples, and overlapping strata cannot be added together. No pooled binomial
confidence interval or universal-equivalence claim is produced. Absence of
differences in this cohort would improve coverage substantially but would not
prove equivalence for all possible variants or resolve ENIGMA's intended setting.

Offline tests cover cohort/reference boundaries, array alignment, analysis,
atomic checkpointing/crash recovery, resume identity and error denominators:

```powershell
.\venv\Scripts\python.exe -m pytest tests/test_spliceai_distance_experiment.py tests/test_spliceai_distance_cases.py tests/test_spliceai_distance_store.py tests/test_spliceai_distance_scan.py -q
```

The original six-SNV report remains unchanged. Its historical script checksum
refers to the archived pilot script in `data/spliceai/build/distance-20260915/`,
not the subsequently extended runner. Neither experiment changes production
SpliceAI settings or classification criteria.

Implementation verification on 16 September 2026: 58 offline tests passed. A
bounded real-model run completed the two fixture SNVs plus a deletion and an
insertion on each strand, then paused at six cases. A separate resumed process
added two cases, yielding eight successful cases and no errors. The outputs are
under `data/spliceai/build/distance-scaled-smoke-20260916/`. This verifies the
runner and resume path, **not** completion of the 4,047-case experiment.

## Full stratified results: 2026-09-17

The complete stratified run finished all **4,047 alleles** with no errors. It
covered 1,315 positions across all 98 canonical BRCA1 and BRCA2 exon and intron
regions: 3,919 SNVs, 64 deletions and 64 insertions. The D=4999, D=5000 and
D=10000 comparison produced no decision-band changes and no maximum-score
differences above the prespecified 1e-6 numerical tolerance.

A separate extension evaluated D=50, D=500 and the shared D=4999 anchor over the
same cohort. Every D=4999 result reproduced the saved base run.

- D=50 versus D=500 changed the SpliceAI decision band for 75 alleles at 52
  distinct loci. There were 44 transitions from 0.10 or lower into the
  intermediate band, 11 from 0.10 or lower to 0.20 or higher, and 20 from the
  intermediate band to 0.20 or higher.
- D=500 versus D=4999 changed 1,106 maximum scores by more than 1e-6, with a
  largest absolute difference of 0.045821, but did not change the decision band
  for any tested allele.
- D=4999, D=5000 and D=10000 had no decision-band differences. The largest
  observed maximum-score difference among these settings was 7.75e-7.

These results show that D=50 is not sufficient for this validation cohort. The
agreement between D=500 and D=4999 is empirical and does not prove equivalence
for every possible variant. Random and targeted cases are mixed, strata overlap,
and three alternate alleles at one locus are correlated. Counts from this design
must not be interpreted as population prevalence estimates.

No production setting or classification rule was changed. ARIANE remains pinned
to D=10000 in `data/spliceai/enigma_v1_2_spliceai_profile.json`.

The formatted report with figures is available as
[`output/pdf/spliceai_distance_experiment_report.pdf`](../../../output/pdf/spliceai_distance_experiment_report.pdf).
