# De novo PD-L1 miniprotein binder design on a single 8 GB GPU

A staged screening pipeline for de novo PD-L1 binder design, run end-to-end on one
**RTX 3070 Ti (8 GB VRAM)** inside WSL2/Ubuntu. It reproduces the multi-predictor
consensus methodology of Anthropic's de novo binder-design study, but **stages it by
compute cost** so it fits in 8 GB: a cheap local co-fold triages a large candidate pool,
and accurate cloud co-folders are spent only on the survivors.

> **Interactive report:** open [`report/index.html`](report/index.html) in a browser —
> full write-up with the screening funnel, result tables, novelty analysis, and
> **interactive 3D predicted complexes** (3Dmol.js) for the top-30 designs.

## Pipeline

```
RFdiffusion            52 backbones      (5C3T / 4ZQK / 4Z18 crystal seeds + de novo scaffolds,
                                          PD-1 epitope hotspots Y56 / R113 / M115 / Y123 fixed)
   │  ProteinMPNN + SolubleMPNN (T=0.2)
   ▼
MPNN sequences         4,905             (2,716 Soluble / 2,189 Protein)
   │  per-backbone cap
   ▼
Candidate pool         2,481             (assigned cand IDs, monomer-folded)
   │  localColabFold monomer foldability, pLDDT ≥ 70
   ▼
Foldable               2,162             (median pLDDT 91.7)
   │  dedup 90% (removed 0) + diversity down-sample 70%
   ▼
Complex co-fold        1,050             (localColabFold AF2-multimer, recycle 1)
   │  Stage-1 ipSAE gate → 242 · Stage-2 + sc_DockQ → 157
   │  liability + redundancy clustering
   ▼
Curated               50
   │  3 cloud co-folders (Boltz-2 + ESMFold2-Full + Protenix)
   │  pose-PASS = all three sc_DockQ ≥ 0.23
   ▼
Pose-agreed           43
   │  composite z = 4·z(ipSAE) + 1·z(sc_DockQ)
   ▼
Reported top          30
```

## Key results

- **43** designs pass all three cloud co-folders; ranked by the dataset's 4:1
  ipSAE : sc_DockQ z-score.
- **Top design `cand00169`**: composite z **+6.04**, ipSAE 0.908, sc_DockQ 0.893.
- **Novelty:** every one of the 43 (and all of the top-30) is **< 45 % identical**
  (max 44.6 %, global BLOSUM62) to any of the 90 reference PD-L1 designs — new
  solutions to the PD-1 epitope, not rediscoveries. For 13 of the top-30 the nearest
  reference neighbour is an experimentally *confirmed* binder.

See `report/index.html` for the full analysis and figures.

## Metrics — computed exactly as in the reference study

- **ipSAE_min** — Dunbrack's `scripts/ipsae.py` (`pae_cutoff 10, dist 10`), min of the two
  `asym` chain directions per structure. Pure-PAE, so invariant to relaxation.
- **sc_DockQ** — DockQ (v2.1.3) of the predicted pose vs the *designed* backbone
  (`scripts/score_iface.py`).
- **ΔΔG / CMS** — PyRosetta `InterfaceAnalyzerMover` after a **coordinate-constrained
  full relax** (`scripts/bench_ddg.py`). See the lesson in `results/LESSONS.md`:
  side-chain-only relax on recycle-1 structures produces false "clash blow-ups".

## Running it

One entry point drives (or documents) every stage in order:

```bash
./run_pipeline.sh list      # list stages
./run_pipeline.sh check      # environment / prereq check
./run_pipeline.sh 07_rank    # run one stage (local stages run for real)
./run_pipeline.sh all        # run local stages, document the GPU/cloud ones
```

See **[PIPELINE.md](PIPELINE.md)** for the full walk-through: the exact command and
environment-switch for every stage, the ColabFold cache/batching tips, and the
field-notes/lessons distributed per stage.

## Repository layout

```
run_pipeline.sh     single entry point / orchestrator
PIPELINE.md         per-stage commands, env switches, caching + batching tips, lessons
inputs/rfdiffusion/ trimmed PD-L1 target structures fed to RFdiffusion (stage 00)
scripts/            pipeline + scoring code
  00_rfdiffusion.sh   stage 00 - the 4 RFdiffusion arms (exact contigs + hotspots)
  cofold/             stage 06 clients: boltz2_run.py, esmfold2_full.py, protenix_run.py
  ipsae.py            ipSAE (vendored from the Dunbrack lab — see header/attribution)
  score_iface.py      sc_DockQ + BSA vs designed backbone
  bench_ddg.py        PyRosetta ΔΔG / CMS with constrained full relax
  rosetta_iface.py    InterfaceAnalyzer helpers
  eval_colabfold.py   localColabFold evaluation
  recompute_pdl1.py   metric recomputation vs the reference dataset
  run_complex.sh      localColabFold complex-fold driver (recycle 1, reused MSA, resumable)
  run_scoring.sh      scoring driver
  analysis/
    rank_designs.py     stage 07 - reproduces final_ranking.csv from the 3 co-fold tables
    align_top30.py      superpose the top-30 complexes onto a common PD-L1 frame
    novelty_and_extract.py / slim_report_data.py / build_report.py   report generation
    plot_zscore_dist.py z-score distribution figure
results/            small result tables (CSV/TXT) for every stage of the funnel
report/index.html   the interactive HTML report (self-contained; 3Dmol.js from CDN)
```

## Reproducing

External tools (not vendored — install separately):

| Stage | Tool |
|---|---|
| Backbone generation | [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) |
| Sequence design | [ProteinMPNN / SolubleMPNN](https://github.com/dauparas/ProteinMPNN) |
| Local co-folding | [localColabFold](https://github.com/YoshitakaMo/localcolabfold) |
| Cloud co-folders | Boltz-2 (`api.boltz.bio`), ESMFold2-Full (Biohub Forge), Protenix (JapanFold) |
| Energetics | [PyRosetta](https://www.pyrosetta.org/) |
| Pose metric | [DockQ](https://github.com/bjornwallner/DockQ) |

Python analysis env: `pip install -r requirements.txt`.
**numpy note:** DockQ 2.1.3 needs `numpy<2`; the ESM SDK needs `numpy 2.x` — pin
`numpy==1.26.4` before scoring, `numpy==2.5.2` before ESM folding.

## Data availability

The reference PD-L1 designs and their measured affinities come from the public
**[Anthropic/claude-protein-binder-design](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design)**
dataset. Those files are **not redistributed here** — download them from HuggingFace.
Bulk intermediates (ColabFold outputs, RFdiffusion backbones, the Python venv) are also
excluded; only the small result tables and the code are tracked.

## Acknowledgements

- Reference methodology and dataset: Anthropic de novo binder-design study.
- ipSAE: Dunbrack lab. DockQ: Wallner lab. Structure viewer: 3Dmol.js.
