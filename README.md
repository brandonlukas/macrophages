# Macrophage GRN Workflow

A Snakemake pipeline to build and interrogate macrophage / dendritic cell transcriptional networks. It combines public ChIP-Atlas binding targets with a curated Zbtb46 list, prunes the network to the supplied cells, infers activity with BITFAM, and runs trajectory-aware analyses and perturbation simulations.

## What it does
- Download TF target sets from ChIP-Atlas (mouse mm10) and merge with curated Zbtb46 targets, then prune the network to expressed genes.
- Infer TF activity with BITFAM (multiple seeds aggregated) using the provided single-cell object.
- Infer pseudotime trees with Lamian, export branches, and run time/variable differential expression on RNA and Z-scores.
- Prepare Seurat objects and compute markers / differential expression for RNA and activity assays.
- Convert data to CellOracle format, build gradients, and simulate knockdown and overexpression perturbations; aggregate perturbation scores and Markov chain transitions.
- Optional PGD-CellOracle branch to repeat simulations on PGD-transformed embeddings.

## Inputs & configuration
- Primary input: `config/config.yaml` expects a `cells` RDS (`inputs/cells.rds` by default) plus resource paths for atlas factors and the Zbtb46 table.
- Tweak parameters (e.g., BITFAM whitelist regex, CellOracle grid size/mass) inside `config/config.yaml`.

## How to run
1. Install Snakemake and ensure the required R/Python environments for the referenced scripts are available (BITFAM, Lamian, Seurat, CellOracle, PGD-CellOracle).
2. Update paths in `config/config.yaml` if your inputs differ.
3. Launch the workflow from the repo root:
   ```bash
   snakemake -s workflow/Snakefile -j 8
   ```

## Key outputs
- `results/BITFAM/network/` and `results/BITFAM/results/`: network files and activity estimates.
- `results/lamian/`: pseudotime, uncertainty, and differential tests.
- `results/seurat/`: markers and DE tables.
- `results/celloracle*/`: CellOracle objects, perturbation scores, and transition summaries for baseline and overexpression scenarios.
- `results/pgd_celloracle*/`: PGD-based variants of the CellOracle outputs.
