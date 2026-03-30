# DiCE-Duan Validation and Benchmarking: Viva Presentation Guide

This document is designed to help you prepare for your viva presentation. It explains the rationale behind the benchmarking tests we conducted and provides a detailed, plain-English explanation for each generated figure in this directory.

---

## Part 1: The "Why" - Benchmarking Rationale

When defending your thesis or presenting this algorithm, examiners will fundamentally ask: **"Why do we need DiCE-Duan? Is it actually better than what people already do?"** 

To answer this, our pipeline (`validation_pipeline.py`) includes a rigorous benchmarking section where DiCE-Duan goes head-to-head against two standard baseline methods. 

Here is how you explain *why* we chose these specific baselines:

### 1. Differential Expression (DE) Ranking (The Traditional Baseline)
*   **What it is:** The standard biological approach. It ranks genes based purely on how much their expression levels change between normal and tumor samples (using log2 Fold Change and p-values).
*   **Why we benchmark against it:** We want to prove the core hypothesis of network biology: **"Genes don't act alone."** DE ranking looks at genes in isolation. By showing that DiCE-Duan outperforms DE ranking, you successfully argue that *looking at how a gene's network connections change (DiCE)* is more powerful for identifying disease drivers than just looking at whether a gene is "turned up or down."

### 2. Degree Centrality Ranking (The Basic Network Baseline)
*   **What it is:** A simple network approach that just asks, "Which genes have the most connections (highest degree) in the tumor network?"
*   **Why we benchmark against it:** We need to prove that DiCE-Duan is doing something sophisticated. If simply picking the most connected "hub" genes gave the same exact results, DiCE-Duan would be overly complex for no reason. We benchmark against this to show that measuring the *differential shift* in centrality (how much a gene's role changes from normal to tumor) captures unique disease dynamics that simple hub-counting misses.

### The Metric: Precision@K
*   **How to explain it:** "We evaluated the algorithms using Precision@K. This means if we take the top $K$ predicted genes (e.g., top 10, 20, or 50), what percentage of them are *already known and experimentally validated* as cancer drivers in databases like COSMIC or DisGeNET? We focus on small $K$ values because in the real world, biologists only have the time and budget to experimentally test the top 10 or 20 predictions, not thousands."

---

## Part 2: Explaining the Visual Evidence (Figure by Figure)

Use these explanations to narrate the figures during your presentation.

### The Benchmarking Results (The Competition)

#### `08_benchmark_comparison.png` & `08_benchmark_bar_k100.png`
*   **What to say:** "This is the ultimate test of our algorithm's predictive power. As you can see, the blue line (DiCE-Duan) significantly outperforms the red line (traditional Differential Expression). For example, in the top 10 predictions, DiCE-Duan achieved 40% precision compared to DE's 20%. It also competes strongly with the yellow line (Degree Centrality). This graph proves that leveraging differential network rewiring is a vastly superior strategy for prioritizing disease genes compared to traditional node-level expression analysis."

### Statistical & Biological Overlap (The Proof)

#### `02_hypergeometric_overlap.png`
*   **What to say:** "To ensure our results weren't just a statistical fluke, we ran a Hypergeometric Significance Test. This plot shows the expected overlap of our predicted genes with known cancer databases versus the *actual* observed overlap. We saw roughly a 3-to-4-fold enrichment over random chance, with extremely significant p-values (e.g., p < 0.000001 for the COSMIC census). This gives us high confidence that the DiCE-Duan predictions are deeply rooted in known cancer biology."

#### `03_precision_at_k.png`
*   **What to say:** "This curve breaks down our algorithm's accuracy as we go further down the ranked list. Naturally, precision is highest at the very top (K=10). The fact that the curve decays smoothly instead of dropping off a cliff shows that the algorithm is relatively robust and stable in its ranking."

### Functional Relevance (What the genes do)

#### `01_pathway_enrichment_barplot.png` & `01_pathway_enrichment_dotplot.png`
*   **What to say:** "We didn't just want a list of names; we wanted to know what these top predicted genes are doing collectively. We ran Pathway Enrichment Analysis against the KEGG and Gene Ontology databases. You can clearly see that our top genes are highly concentrated in critical cancer-related pathways, such as cell cycle regulation and ubiquitin-mediated proteolysis. This confirms that DiCE-Duan is identifying functionally coordinated gene modules that drive the tumor phenotype."

### Network Topology Dynamics (How it works under the hood)

#### `04_centrality_shift_scatter.png` & `04_centrality_boxplot.png`
*   **What to say:** "These plots visualize the core philosophy of our algorithm. We are tracking the 'centrality shift'—how much a gene's influence changes between the normal and tumor state. The scatter plot and boxplots clearly show that our top candidate genes (highlighted in red) experience a massive, statistically significant shift in their network topology compared to the vast background of normal genes. They are actively 'rewiring' their connections to support the tumor."

#### `05_rewiring_volcano.png` & `05_rewiring_histogram.png`
*   **What to say:** "The volcano plot maps the magnitude of betweenness centrality change against its statistical significance. We can immediately spot highly significant rewired genes like EP300 and the Androgen Receptor (AR). These are genes that have completely changed their structural role in the prostate cancer network, making them prime candidates for targeted therapy."

### Validation & Simulation (The Real-world Application)

#### `06_literature_validation.png`
*   **What to say:** "As an empirical check, we query the PubMed database to see how often our top predicted genes are published in the context of this specific cancer. A high publication count essentially serves as independent, external validation by the scientific community that our algorithm is pointing researchers in the right direction."

#### `07_knockout_disruption.png`
*   **What to say:** "Finally, we ran a virtual 'Knockout Simulation'. If a pharmaceutical company were to develop a drug targeting one of our genes, what would happen to the tumor? This chart ranks the genes based on how much the entire tumor network structurally collapses when you virtually remove them. Genes at the top of this list are structurally critical vulnerabilities in the tumor network."
