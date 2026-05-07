# Diagnostic Dossier

- **Commit**: `fe3000e54752b56b22e625729ce0bc6698b0d207`
- **Generated**: 2026-05-07T11:46:31.465711+00:00
- **N probes**: 9

---

## Stage 1-raw-inputs — `raw_inputs`

- _duration_: 0.907s
- **expression_normal**:
  - **path**: data/raw/expression/expression_normal.txt
  - **n_genes**: 10868
  - **n_samples**: 52
  - **nan_cells**: 0
  - **duplicate_genes**: 0
  - **id_type**: symbol
- **expression_tumor**:
  - **path**: data/raw/expression/expression_tumor.txt
  - **n_genes**: 10868
  - **n_samples**: 498
  - **nan_cells**: 0
  - **duplicate_genes**: 0
  - **id_type**: symbol
- **kegg**:
  - **path**: data/raw/kegg
  - **n_xml_files**: 369

## Stage 10-dnf — `dnf`

- _duration_: 0.000s

## Stage 2-parsers — `parsers`

- _duration_: 1.347s
- **kegg**:
  - **path**: data/processed/networks/kegg_interactions.txt
  - **n_edges**: 6222
  - **duplicate_edges**: 0
  - **unique_genes**: 2530
  - **self_loops**: 0
- **biogrid**:
  - **path**: data/processed/networks/biogrid_interactions.txt
  - **n_edges**: 983026
  - **duplicate_edges**: 0
  - **unique_genes**: 20361
  - **self_loops**: 0

## Stage 3-merge — `merge_networks`

- _duration_: 3.879s
- **path**: data/processed/networks/merged_network.txt
- **n_edges**: 1856918
- **n_directed**: 6222
- **n_undirected**: 1850696
- **unique_genes**: 20487
- **bidirectional_pairs**: 171
- **self_loops**: 0

## Stage 4-filters — `filters`

- _duration_: 0.032s
- **dea**:
  - **n_genes_total**: 10868
  - **n_significant**: 1269
  - **fdr_p10**: 5.85049e-12
  - **fdr_p50**: 0.00192101
  - **fdr_p90**: 0.629331
  - **fdr_le_0_05_pct**: 0.655226
  - **abs_lfc_gt_1_pct**: 0.116765
  - **n_genes_listed**: 1269
- **ig**:
  - **n_genes**: 4129
- **intersection**:
  - **n_genes**: 1165
  - **n_dea_only**: 104
  - **n_ig_only**: 2964
  - **jaccard_dea_ig**: 0.275219

## Stage 5-weighting — `weighting`

- _duration_: 4.879s
- **n_input_edges**: 1856918
- **n_output_edges**: 490436
- **n_dropped**: 1366482
- **edges_unresolved_in_expression**: 1366482
- **gene_symbol_resolution_rate**: 0.448577
- **epsilon_clamped_edges**: 2
- **weight_truncation_collisions**: 285839
- **weight_int_eq_1_count**: 2
- **weight_distribution**:
  - **min**: 1e-06
  - **max**: 1
  - **p10**: 0.471242
  - **p50**: 0.783803
  - **p90**: 0.962337

## Stage 6-centrality — `centrality`

- _duration_: 0.010s
- **omp_determinism_max_relative_drift**: 0
- **sampled_vs_full_subgraph_n**: 1000
- **sampled_vs_full_spearman_rho**: 0.9315
- **sampled_vs_full_spearman_p**: 0

## Stage 7-differential — `differential`

- _duration_: 0.017s
- **n_genes**: 8960
- **bet_eig_rank_correlation**: 0.621658
- **top12_eig_dropped_n**: 1128
- **top12_eig_dropped_genes**: [TNF, PSME3, PRKACA, CTNNB1, HIF1A, NTRK1, SMURF1, PRNP, HSPA5, HTT...]
- **top12_eig_dropped_in_cgc**: 121
- **cgc_size**: 685

## Stage 8-validation — `validation`

**ERROR**: `ParserError: Error tokenizing data. C error: Expected 1 fields in line 8, saw 4
`


---

## Bug Registry

_Populate from per-stage findings; sort by impact × confidence × blast_radius._

| ID | Stage | Symptom | Evidence | Confidence | Blast radius | Hypothesized root cause |
|----|-------|---------|----------|------------|--------------|-------------------------|
