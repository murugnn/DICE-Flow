# DiCE-Duan Python Codebase — Complete Reference

This document explains every Python file in the DiCE-Duan project, organized by pipeline phase.

---

## Pipeline Scripts (`src/python/`)

These are the core backend modules that execute sequentially via `run_pipeline.sh`.

### Phase 1: Network Construction

#### `parse_kegg.py`
- **Purpose:** Extracts directed gene–gene interactions from KEGG pathway XML (KGML) files.
- **How it works:** Uses Python's `xml.etree.ElementTree` to parse each `.xml` file. It builds an `entry_map` mapping KEGG entry IDs to human gene symbols (extracted from the `<graphics>` tag). It then iterates over `<relation>` elements (types `PPrel`, `GErel`, `PCrel`) to extract directed edges like `TP53 → MDM2`. Self-loops are filtered out.
- **Key Class:** `KEGGParser` — initialized with a directory of KEGG XML files.
- **Output:** A TSV file of `(source, target, direction)` tuples.

#### `parse_biogrid.py`
- **Purpose:** Parses the massive BioGRID TAB2 experimental interaction database, filtering for human (Homo sapiens, Taxonomy ID 9606) protein–protein interactions.
- **How it works:** Reads the BioGRID file in **chunks of 100,000 rows** to avoid memory overload. For each chunk, it filters by organism ID on both interactors, extracts gene symbols, removes self-loops and duplicates.
- **Key Class:** `BioGRIDParser` — handles format detection and chunked I/O.
- **Output:** A TSV file containing all unique human PPIs with columns `(source, target, interaction, type)`.

#### `merge_networks.py`
- **Purpose:** Combines KEGG (directed) and BioGRID (undirected) interactions into a single master scaffold network.
- **How it works:** Loads both parsed files. KEGG edges are stored as directed. BioGRID edges that are *not* already represented in KEGG are added bidirectionally (A→B and B→A). Supports an optional `--candidates` filter to restrict the network to a pre-filtered gene list.
- **Key Class:** `NetworkMerger` — handles deduplication and source-database tracking.
- **Output:** `merged_network.txt` — the master edge list.

### Phase 2: Expression Integration

#### `pre_filter.py`
- **Purpose:** Implements DiCE Phase II — reduces the gene pool using **Information Gain** (Mutual Information) to remove genes that carry no discriminatory signal between Normal and Tumor states.
- **How it works:** Loads both expression matrices, transposes them (patients as rows), and creates a binary label vector (0 = Normal, 1 = Tumor). Uses `sklearn.feature_selection.mutual_info_classif` to compute each gene's Information Gain score. Genes scoring above the mean IG threshold are kept.
- **Output:** A text file listing one gene symbol per line — the filtered candidate set.

#### `weight_network.py`
- **Purpose:** Maps biological context onto the structural scaffold by calculating **Pearson Correlation** weights for every edge.
- **How it works:** For each edge (gene_A, gene_B), it loads their expression vectors across all patients and computes the Pearson correlation coefficient `r`. The DiCE edge weight formula is applied: `W = 1 − |r|`. Perfectly correlated genes (|r| = 1) get a weight of 0 (shortest distance); independent genes (r = 0) get weight 1. Weights are also scaled to integers (×1000) for the C++ Radix Heap engine. A minimum epsilon of 1e-6 prevents zero-weight edges.
- **Key Class:** `NetworkWeighter` — supports auto-detection of expression file formats (CSV, TSV, GCT).
- **Output:** `network_normal.txt` and `network_tumor.txt` — condition-specific weighted edge lists.

### Phase 3: Scoring

#### `differential_centrality.py`
- **Purpose:** Implements DiCE Phase V — the **Ensemble Ranking** that produces the final disease gene list.
- **How it works:**
  1. Loads the centrality outputs from the C++ engine for Normal and Tumor networks.
  2. Calculates **absolute differences** (|ΔCB| and |ΔCE|) between conditions.
  3. Converts to **percentile ranks** (0.0 to 1.0) using Pandas `.rank(pct=True)`.
  4. Computes the **Rank Product** = `rank_bet × rank_eig` — penalizing genes that only shifted in one metric.
  5. Applies a **noise filter:** genes with below-average betweenness in *both* conditions are excluded (they are peripheral nodes where small shifts are mathematical artifacts).
  6. Adds signed delta columns for downstream visualization.
- **Key Class:** `DifferentialCentralityAnalyzer`.
- **Output:** `dice_genes.txt` (full) and `dice_genes_top500.txt` (top 500).

### Phase 4: Validation & Visualization

#### `validate_candidates.py`
- **Purpose:** A lightweight sanity-check script that computes **raw degree** (number of connections) for the top 20 predicted genes directly from the weighted network files.
- **How it works:** Loads the network edge lists, counts how many times each gene appears (as source or target), and prints a table showing Normal Degree vs Tumor Degree with a human-readable status (`Gained Hubness`, `Lost Hubness`, `Low Confidence`).
- **Output:** Console output table.

#### `export_results.py`
- **Purpose:** Generates **publication-quality figures** and exports data for **Cytoscape** network visualization.
- **Key Features:**
  - `prepare_cytoscape_network()` — formats weighted networks as Cytoscape-importable CSVs.
  - `prepare_node_attributes()` — exports DiCE scores and regulation direction as Cytoscape node attribute files.
  - `plot_dice_rankings()` — horizontal bar chart of top-N genes by ensemble score.
  - `plot_centrality_scatter()` — log-scale scatter plot of Normal vs Tumor betweenness with top-10 gene labels.
- **Key Class:** `NetworkVisualizer`.

#### `generate_dummy_data.py`
- **Purpose:** A testing utility that generates **synthetic** KEGG XML files, BioGRID interaction tables, and expression matrices. Used by the `--test` flag in `run_pipeline.sh` to verify pipeline connectivity without waiting for real TCGA data.
- **Output:** A complete `data/dummy/` directory mimicking the real data structure.

### Web Interface

#### `app.py`
- **Purpose:** A **Streamlit** interactive dashboard for exploring DiCE results visually.
- **Key Features:**
  - **Feature 1 — Path Switch:** Users input a source and target gene. The app computes the shortest signaling path in both Normal and Tumor networks using NetworkX, displays path costs, and renders interactive PyVis network graphs.
  - **Feature 2 — Rewiring Proof:** Users inspect a candidate gene's "ego network." The app compares its average edge weight (connection strength) between Normal and Tumor, proving whether the gene gained or lost connectivity.
- **Dependencies:** `streamlit`, `networkx`, `pyvis`.

#### `server.py`
- **Purpose:** A **FastAPI** REST API backend exposing the DiCE results as programmatic endpoints.
- **Key Endpoints:**
  - `GET /api/network/initial` — Returns the top-N DiCE genes as a D3-compatible JSON graph.
  - `GET /api/path?source=X&target=Y` — Computes shortest paths in both networks and returns cost difference.
  - `GET /api/rewiring?gene=X` — Returns the ego network for a gene, classifying each edge as `maintained`, `lost`, or `gained`.
  - `GET /api/knockout` — Returns virtual knockout disruption scores.
- **Dependencies:** `fastapi`, `uvicorn`, `networkx`.

---

## Root-Level Scripts

#### `validation_pipeline.py` (930 lines)
- **Purpose:** The **monolithic validation and benchmarking framework** — the largest and most critical Python file. Runs 8 independent scientific analyses to mathematically prove the output quality.
- **Analyses:**
  1. **Pathway Enrichment** — Uses `gseapy` (Enrichr API) against KEGG and GO databases. Generates barplots and dot plots of significant terms.
  2. **Hypergeometric Overlap** — Tests whether the overlap between DiCE's top 100 predictions and known cancer gene databases (COSMIC, DisGeNET) is statistically significant using `scipy.stats.hypergeom`.
  3. **Precision@K** — Evaluates how many of the top K predictions (K=10,20,50,100,200,500) are true positives.
  4. **Network Topology Validation** — Computes centrality shift scatter plots and Mann-Whitney U boxplots proving DiCE candidates are statistical outliers.
  5. **Differential Rewiring** — Generates volcano plots mapping ΔBetweenness vs z-score significance.
  6. **Literature Validation** — Queries NCBI PubMed via `Biopython.Entrez` for publication counts of each top gene in the context of prostate cancer.
  7. **Virtual Knockout Simulation** — Visualizes the C++ engine's knockout disruption scores.
  8. **Benchmarking** — Computes Differential Expression and Degree Centrality rankings as baselines, then plots head-to-head Precision@K comparisons.
- **Output:** Publication-quality PNG figures in `data/results/validation/figures/` and CSV statistical tables in `data/results/validation/tables/`.

#### `survival_validation.py` (292 lines)
- **Purpose:** Clinical survival validation module that bridges computational predictions with real patient outcomes.
- **How it works:**
  1. Loads patient expression and clinical metadata (days_to_death, event).
  2. For each top DiCE gene, stratifies patients into High/Low expression groups (median split).
  3. Performs **Kaplan-Meier** survival curve estimation with confidence intervals.
  4. Runs a **Log-Rank Test** to determine if survival differs significantly between groups.
  5. Fits a **Cox Proportional Hazards** model to compute Hazard Ratios and p-values.
  6. Supports TCGA-style barcode matching (truncation to 12 chars).
- **Dependencies:** `lifelines`.
- **Output:** Kaplan-Meier curve PNGs and `survival_results.csv` with Hazard Ratios and p-values.

#### `post_analysis.py` (224 lines)
- **Purpose:** Post-hoc visualization and gene list export utilities.
- **Key Features:**
  - `DiCEVisualizer` class: Generates top-gene bar charts, centrality change scatter plots, score distribution histograms, and centrality heatmaps.
  - `EnrichmentAnalyzer` class: Exports gene lists formatted for external enrichment tools (Enrichr, DAVID, g:Profiler).
- **Output:** PNG figures and plain-text gene lists for web-based enrichment uploads.
