# De novo PD-L1 miniprotein binder design on a single 8 GB GPU

A staged, compute-budgeted screening pipeline for de novo PD-L1 binder design, run
end-to-end on one **RTX 3070 Ti (8 GB VRAM)** inside WSL2/Ubuntu. The research question
is deliberately narrow and engineering-flavoured:

> **Can a cheap local co-fold pass plus a few free cloud co-folders reproduce an expensive
> multi-predictor consensus for binder triage on hardware that cannot run the full stack?**

> ⚠️ **Computational proof-of-concept — not experimentally validated.** Every "binder" here
> is a *design and an in-silico prediction*. Nothing in this repository has been expressed,
> purified, or assayed. No claim is made that any design actually binds PD-L1. See
> [Limitations](#limitations).

> **Interactive report:** open [`report/index.html`](report/index.html) in a browser —
> the full write-up with the screening funnel, validation figures, result tables, a novelty
> analysis, and **interactive 3D predicted complexes** (3Dmol.js) for the top-30 designs.

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
   │  Stage-1 (ipSAE + sc_DockQ) → 242 · Stage-2 ΔΔG/CMS physics, PyRosetta → 157
   │  (PyRosetta removes 85 weak/repulsive interfaces: ΔΔG ≤ −40 REU & CMS ≥ 360 Å²)
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

## Does the cheap local pass actually work? (validation)

The pipeline only makes sense if a single recycle-1 localColabFold prediction is a valid
triage filter. Two checks (see `report/index.html` §02 and `docs/figures/`):

- **Metric code is correct.** Re-computing ipSAE and sc_DockQ on the reference dataset's *own*
  structures reproduces its published values (r ≈ 1.0) — so the same code applied to our
  predictions measures the same thing.
- **localColabFold tracks the consensus.** On 7 reference PD-L1 designs, our independent
  localColabFold ipSAE_min / ipTM agree with the dataset's 10-model consensus
  (Pearson r = 0.968 / 0.981, n = 7), and the wetlab non-binders that the consensus scores
  low are scored low by localColabFold too. A single cheap pass already separates the known
  binders from the known non-binders in this small set.

## Key results (in silico)

- **43** of the 50 curated designs pass **all three** cloud co-folders at `sc_DockQ ≥ 0.23`;
  they are ranked by the dataset's 4:1 `ipSAE : sc_DockQ` z-score.
- **Top-ranked design `cand00169`**: composite z ≈ +6.0 (co-fold ipSAE ≈ 0.91, sc_DockQ ≈ 0.89).
- **Sequence novelty:** each of the 43 (and all of the top-30) is **< 45 % identical**
  (max 44.6 %, global BLOSUM62) to any of the 90 reference PD-L1 designs — i.e. these are
  *new sequences*, not copies of the reference set. This is a **sequence-level** statement
  only; it is **not** evidence that the designs bind.
- For 13 of the top-30, the *nearest* (still < 45 % identical) reference sequence happens to
  be an experimentally confirmed binder. This is an observation about *where in sequence
  space the designs sit*, **not** a functional result for our designs.

## Metrics — computed as in the reference study

- **ipSAE_min** — Dunbrack's `scripts/ipsae.py` (`pae_cutoff 10, dist 10`), min of the two
  `asym` chain directions per structure. Pure-PAE, so invariant to relaxation.
- **sc_DockQ** — DockQ (v2.1.3) of the predicted pose vs the *designed* backbone
  (`scripts/score_iface.py`).
- **ΔΔG / CMS** — PyRosetta `InterfaceAnalyzerMover` after a **coordinate-constrained full
  relax** (`scripts/bench_ddg.py`); see `results/LESSONS.md` (side-chain-only relax on
  recycle-1 structures produces spurious "clash blow-ups").

## Running it

```bash
export PDL1_PROJECT=/path/to/working/data   # dir holding structures / scoring CSVs / dataset
export PDL1_OUT=./build                      # where generated report data + figures go
./run_pipeline.sh list       # list stages
./run_pipeline.sh check       # environment / prereq check
./run_pipeline.sh 07_rank     # run one stage (local stages run for real on the committed CSVs)
```

**Configuration.** No paths are hard-coded: every script reads `PDL1_PROJECT` (working data
directory) and `PDL1_OUT` (output directory) from the environment, with a `~`-relative default.
`run_pipeline.sh` and cloud clients also take explicit CLI flags. See **[PIPELINE.md](PIPELINE.md)**
for the exact command and environment-switch for every stage (including the numpy 1.26.4 ↔ 2.5.2
toggle and the ColabFold cache/batching tips).

The pure-Python stages run from a clean clone against the committed tables, e.g.:
```bash
python scripts/analysis/rank_designs.py --scoring results --verify   # reproduces final_ranking.csv
```

## Repository layout

```
run_pipeline.sh     single entry point / orchestrator (env-driven, no hard-coded paths)
PIPELINE.md         per-stage commands, env switches, caching + batching tips, lessons
inputs/rfdiffusion/ trimmed PD-L1 target structures fed to RFdiffusion (stage 00)
scripts/            pipeline + scoring code
  00_rfdiffusion.sh   stage 00 — the 4 RFdiffusion arms (exact contigs + hotspots)
  cofold/             stage 06 clients: boltz2_run.py, esmfold2_full.py, protenix_run.py
  ipsae.py            ipSAE (vendored from the Dunbrack lab — see header/attribution)
  score_iface.py      sc_DockQ + BSA vs designed backbone
  bench_ddg.py        PyRosetta ΔΔG / CMS with constrained full relax
  run_complex.sh / run_scoring.sh   localColabFold complex-fold + scoring drivers
  analysis/           filter_stages.py (Stage-1/Stage-2 gates), rank_designs.py, align_top30.py,
                      make_valfig.py, and the report generator
results/            small result tables (CSV/TXT) for every stage of the funnel
docs/figures/       validation figures
report/index.html   the interactive HTML report (self-contained; 3Dmol.js from CDN)
```

## Reproducing (external tools)

Not vendored — install separately:

| Stage | Tool |
|---|---|
| Backbone generation | [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) |
| Sequence design | [ProteinMPNN / SolubleMPNN](https://github.com/dauparas/ProteinMPNN) |
| Local co-folding | [localColabFold](https://github.com/YoshitakaMo/localcolabfold) |
| Cloud co-folders | Boltz-2 (`api.boltz.bio`), ESMFold2-Full (Biohub Forge), Protenix (JapanFold) |
| Energetics | [PyRosetta](https://www.pyrosetta.org/) |
| Pose metric | [DockQ](https://github.com/bjornwallner/DockQ) |

Python analysis env: `pip install -r requirements.txt`. **numpy note:** DockQ 2.1.3 needs
`numpy<2`; the ESM SDK needs `numpy 2.x` — pin `numpy==1.26.4` for scoring, `numpy==2.5.2`
for ESM folding.

## Limitations

- **In silico only.** No experimental expression, purification, or binding assay. All results
  are model predictions; treat "pass / top-ranked" as *computational triage*, not validated hits.
- **Reproduction, not new methodology.** The selection metrics and thresholds
  (ipSAE_min, `sc_DockQ ≥ 0.23`, 4:1 z-weighting, top-30/target) are **adopted from the
  Anthropic study**, not derived or independently calibrated here.
- **Consensus is a subset.** The 3-way consensus uses Boltz-2 + ESMFold2-Full + Protenix —
  the co-folders that were free/accessible on this hardware — not the dataset's full 10-model panel.
- **Small validation set.** The consensus-agreement correlations are on n = 7 designs.
- **External novelty skipped.** Novelty is measured only against the reference PD-L1 set, not
  against UniRef90 (the intended MMseqs2/UniRef90 screen was blocked by unreachable servers and
  bandwidth limits). "Novel" here means *distinct from the reference designs*, not *absent from
  known proteins*.
- **Boltz-2 is stochastic** (`num_samples=1`), so its per-design metrics vary run to run.

## What's mine vs. adopted

- **This project (my work):** the compute-staged pipeline design and orchestration; running
  RFdiffusion + ProteinMPNN/SolubleMPNN + localColabFold locally on an 8 GB GPU to generate,
  fold, and filter the pool; the scoring/ranking implementation; the three cloud co-folder API
  clients; the validation experiments; and the interactive report.
- **Adopted / third-party:** the overall selection methodology, metrics, and thresholds and the
  reference designs + measured affinities come from the public
  **[Anthropic/claude-protein-binder-design](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design)**
  study and dataset (not redistributed here — download from HuggingFace). RFdiffusion,
  ProteinMPNN, localColabFold, Boltz-2, ESMFold2, Protenix, PyRosetta, DockQ, and ipSAE are
  third-party tools. This repository was built with substantial AI-assisted coding (Claude Code),
  directed and reviewed by me.

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgements

Reference methodology and dataset: Anthropic de novo binder-design study. ipSAE: Dunbrack lab.
DockQ: Wallner lab. Structure viewer: 3Dmol.js.
