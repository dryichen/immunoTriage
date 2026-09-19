# immunoTriage

Calibrating immune gene-expression signatures against chance, and predicting how a
signature will behave from the cells that produce its signal.

## What it does

A signature that separates outcomes is not thereby informative: random gene sets of the
same size, scored the same way, often do as well. `immunoTriage` measures every signature
against a random-gene-set null matched on **gene-set size and scoring engine**, scores each
signature with the method of its own publication, and tests whether an immune-specific
result survives confounder controls (proliferation, hypoxia, chromosomal instability).

In the accompanying analysis of 52 published signatures, the inflamed/interferon family
gained +1.54 null standard deviations for checkpoint response against +0.27 for survival
(permutation P = 0.004): the same signatures that predict response carry little prognostic
information without the treatment.

The package also asks which cells produce a signature's signal, using any single-cell
reference atlas with cell-type labels, and predicts from that alone whether the signature
is worth testing in a checkpoint cohort:

- trained on 50 signatures and a gastric atlas: held-out Spearman rho = 0.68, AUC 0.92 for
  identifying a signature that works;
- predictions locked on 2026-09-14 and tested, under a registration deposited before
  scoring, on a checkpoint-treated cohort never used in development (metastatic urothelial
  carcinoma, GSE176307): rho = +0.38, one-sided P = 0.003; 15 of the 20 signatures called
  for testing worked there, against 7 of 19 called not worth a cohort;
- rebuilt from a lung adenocarcinoma atlas (GSE131907), the rule predicts almost as well
  (held-out rho = 0.65), although the registered permutation test of atlas independence
  narrowly failed (P = 0.072). Treat the shipped model as calibrated on the gastric atlas.

## Install

```r
# install.packages("remotes")
remotes::install_github("dryichen/immunoTriage")
```

Optional engines: `GSVA` (Bioconductor) for GSVA and ssGSEA scoring; `ranger` and
`xgboost` for the ceiling estimate.

## Before you run a cohort: attribute and triage

```r
library(immunoTriage)

# atlas_expr: genes x cells, log-normalized; atlas_type: one label per cell
att <- attribute_signature(atlas_expr, atlas_type, my_signatures)
triage_signature(att)          # predicted response effect and a recommendation

# or everything at once, adding the matched-null benchmark if a cohort is supplied
res <- signature_workflow(my_signatures, atlas_expr, atlas_type,
                          expr = bulk_expr, outcome = response)
```

The shipped model expects the ten classes it was trained on (T, NK, B, plasma,
macrophage, DC, fibroblast, endothelium, pericyte, epithelium); pass `cell_map` if your
atlas uses other labels.

## With a cohort: benchmark against the matched null

```r
library(immunoTriage); library(survival)
sigs <- c(gastro_signatures, gastro_controls)

prog <- benchmark_signatures(gc_tcga$expr, Surv(gc_tcga$time, gc_tcga$event), sigs)
resp <- benchmark_signatures(gc_ici$expr, gc_ici$response, sigs)
test_dissociation(prog, resp, gastro_families)
```

On the bundled gastric data the immune family shifts by +0.26 between settings
(P = 0.004) while the confounder controls do not (interaction P = 0.03).

## Functions

| function | purpose |
|---|---|
| `score_signatures()` | score with each signature's published method (GSVA, ssGSEA, weighted mean, z-mean) |
| `benchmark_signatures()` | compare against a random-gene-set null matched on size and engine, with post hoc direction selection |
| `test_dissociation()` | prognostic against predictive shift, with a mandatory confounder family |
| `test_interaction()` | signature × treatment interaction in a randomized trial, with the random-set false-positive rate measured rather than assumed |
| `attribute_signature()` | which cell types produce a signature's signal, in any annotated single-cell atlas |
| `triage_signature()` | predicted checkpoint-response effect and a recommendation, from attribution alone |
| `refine_signature()` | drop the genes that report on the wrong cells |
| `estimate_ceiling()`, `compare_to_ceiling()` | what a full-transcriptome model reaches on the same data, and how much headroom a signature leaves |
| `signature_workflow()` | attribution, triage, repair and (optionally) benchmark in one call |

## Data

| object | contents |
|---|---|
| `gc_ici` | pembrolizumab-treated gastric cancer (Kim et al. 2018, PRJEB25780) |
| `gc_tcga` | TCGA-STAD primary tumors with overall survival |
| `gastro_signatures`, `gastro_controls`, `gastro_families` | example signatures, confounder controls and their families |
| `triage_model` | the fitted ridge model, its two cut-offs and validation statistics |
| `sig_behavior` | measured prognostic and response effects of the 50 training signatures |
| `inst/extdata/` | the locked predictions, and the attribution of every signature in the gastric and lung atlases |

## Design notes

The random null allows direction to be chosen post hoc, the latitude a published
signature receives in practice, which makes it conservative rather than lenient.
`test_dissociation()` requires a confounder family because a shift in the immune family is
uninterpretable unless the confounders stay put. Signatures are correlated, so
significance across many of them should be judged against a permutation that keeps that
correlation rather than a binomial.

## Reproducing the paper

The analysis scripts, the registrations of both external tests and the derived tables are
in the companion repository `dryichen/immunoTriage-paper`.

## Citation

immunoTriage: calibrating immune transcriptomic signatures against chance and resolving their
cellular origin. Manuscript in preparation.

## License

MIT
