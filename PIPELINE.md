# PIPELINE — running the PD-L1 binder design end to end

Single entry point: [`run_pipeline.sh`](run_pipeline.sh).

```bash
./run_pipeline.sh list      # list stages
./run_pipeline.sh check     # environment / prereq check for every stage
./run_pipeline.sh 07_rank   # run one stage
./run_pipeline.sh all       # run local stages, document the GPU/cloud ones
```

> This is **not** a push-button reproduction on a fresh machine — the pipeline needs
> external tools (RFdiffusion, ProteinMPNN, localColabFold, PyRosetta), a GPU, and
> cloud API keys. `run_pipeline.sh` is the orchestrator that defines the exact order,
> runs the locally-runnable stages for real (05 score, 07 rank, 08 report), and for the
> heavy stages prints the exact command + checks prerequisites. Every stage also prints
> the environment-switch command it needs, so you never have to guess.

Field-tested and verified in this environment: `check` (all green), `07_rank`
(reproduces `final_ranking.csv` to |Δ| = 0.0005), `05_score` on cand00169
(ipSAE + `sc_DockQ 0.892`).

---

## The two environments (and how to switch)

No paths are hard-coded. Set two environment variables (both default to a `~`-relative path):

```bash
export PDL1_PROJECT=/path/to/working/data   # structures, scoring CSVs, the reference dataset
export PDL1_OUT=./build                       # generated report data + figures
```

`run_pipeline.sh` reads `PROJECT` (defaulting to `$HOME/protein_designs/pdl1_bench`) and
exports `PDL1_PROJECT` / `PDL1_OUT` so every Python subprocess inherits them; each Python
script also honours these variables directly.

### ROOT env — RFdiffusion + localColabFold (stages 00, 02, 04)
Installed under `/root/protein`, so run them **as root**, GPU required:
```bash
sudo -i
source /root/miniforge3/etc/profile.d/conda.sh
export PATH=/root/protein/localcolabfold/colabfold-conda/bin:$PATH
export JAX_COMPILATION_CACHE_DIR=$HOME/.cache/colabfold_jax   # reuse compiled kernels
unset DISPLAY; export MPLBACKEND=Agg
```

### USER venv — scoring / ranking / analysis (stages 05, 07, 08)
```bash
source $PDL1_PROJECT/.venv/bin/activate
# or call the interpreter directly: $PROJECT/.venv/bin/python
```

### The numpy toggle (important)
The same venv is used for DockQ scoring and for the ESM SDK, but they need
**different numpy majors**:
```bash
pip install -q numpy==1.26.4    # BEFORE stage 05 / 07 (DockQ 2.1.3 needs numpy<2)
pip install -q numpy==2.5.2     # BEFORE the ESMFold2 client in stage 06 (ESM SDK needs numpy 2.x)
```
`run_pipeline.sh` calls `ensure_numpy` to flip this automatically per stage.
> **Tip — the numpy tug-of-war.** PyPI is fast, so the toggle costs seconds; don't try
> to keep one numpy for both — DockQ silently returns empty `sc_dockq` under numpy 2.x.

### Cloud keys — from the environment, never hardcoded
```bash
export BOLTZ_API_KEY=sk_bc_ws_live_...   # api.boltz.bio (a live key; test keys return a MOCK fixture)
export BIOHUB_TOKEN=...                  # ESMFold2-Full via Biohub Forge
# Protenix (JapanFold) needs no key.
```

---

## ColabFold: the two time-savers (stages 02 & 04)

1. **Reusable target MSA.** PD-L1's MSA is computed **once** and stored in
   `$PROJECT/fold_screen/target_msa/` (`uniref.a3m`, `bfd.*.a3m`, `pdb70.m8`). It is then
   *paired* into every binder's `.a3m`, so the expensive MSA search runs once instead of
   1,050 times. Build/refresh it with the MSA step in `target_msa/msa.sh`.
2. **JAX compilation cache.** `export JAX_COMPILATION_CACHE_DIR=$HOME/.cache/colabfold_jax`
   (~84 MB). The first fold JIT-compiles the model; every later fold reuses the cached
   kernels and starts in seconds. Point every ColabFold run at the same dir.

### Batching tips
- **Split the pool into chunks.** The 1,050 complexes run as 4 chunks (`chunk0..3`).
  [`run_complex.sh`](scripts/run_complex.sh) always attacks the lowest incomplete chunk,
  writes a `.done.txt` per finished design, and **auto-restarts ColabFold on SIGKILL** —
  so any crash/restart resumes instead of starting over.
- **Parallel scoring.** Interface scoring (CPU) parallelises with N independent workers
  each taking a static slice `cands[K::N]`, writing one row-file per candidate (resumable).
- **WSL: poll sparingly.** Aggressive `wsl.exe` polling that times out (SIGTERM) can
  destabilise and restart the WSL VM — which **wipes `/tmp`**. Poll long runs every
  20–30 min, keep workers resumable, and never store anything you need in `/tmp`.

---

## Stages

Overarching principle — **stage by cost, not by taste.** 8 GB can't fold everything well:
run one cheap local co-fold on all 1,050, then spend accurate cloud compute only on the
~50 that survive. The expensive models never see a hopeless design.

### 00 · RFdiffusion → 52 backbones  *(ROOT env, GPU)*
Inputs are committed: the trimmed PD-L1 target structures in **`inputs/rfdiffusion/`**
(`pdl1_target.pdb`, `pdl1_5C3T.pdb`, `pdl1_4ZQK.pdb`, `pdl1_4Z18.pdb`). Run all four arms:
```bash
sudo -i && source /root/miniforge3/etc/profile.d/conda.sh && conda activate <rfdiffusion env>
bash scripts/00_rfdiffusion.sh          # loops the 4 arms below; set NUM=/OUT= to taste
```
Exact settings, recovered from the run `.trb` configs (checkpoint `Complex_base_ckpt.pt`, zero noise):

| arm | input_pdb | contig | binder len | hotspot_res |
|---|---|---|---|---|
| pdl1_binder (de novo) | `pdl1_target.pdb` | `A2-110/0 70-80` | 70–80 | `A39,A96,A98,A106` |
| 5C3T_binder | `pdl1_5C3T.pdb` | `A18-132/0 65-85` | 65–85 | `A56,A113,A115,A123` |
| 4ZQK_binder | `pdl1_4ZQK.pdb` | `A18-132/0 65-85` | 65–85 | `A56,A113,A115,A123` |
| 4Z18_binder | `pdl1_4Z18.pdb` | `A18-132/0 65-85` | 65–85 | `A56,A113,A115,A123` |

One invocation used `num_designs` 4–8; the script was looped to build ~50 backbones per template,
then diversity-reduced to 52. The de-novo target is numbered 2–110 and the crystal templates 18–132,
so the two hotspot sets mark the **same epitope** offset by 17.
```bash
# the arm each row expands to:
python /root/protein/RFdiffusion/scripts/run_inference.py \
  inference.input_pdb=inputs/rfdiffusion/pdl1_5C3T.pdb inference.output_prefix=$OUT/5C3T_binder \
  inference.num_designs=8 'contigmap.contigs=[A18-132/0 65-85]' 'ppi.hotspot_res=[A56,A113,A115,A123]' \
  denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
```
> **Tip — seed from a real epitope.** Backbones grown against the PD-1-binding face from
> crystal complexes **5C3T / 4ZQK / 4Z18** plus de novo scaffolds, hotspots
> **Y56 / R113 / M115 / Y123** (crystal numbering) fixed. Anchoring on a known epitope keeps
> the whole downstream funnel pointed at a bindable surface.

### 01 · ProteinMPNN + SolubleMPNN → 4,905 sequences  *(MPNN env, CPU ok)*
Run **both** designers at sampling temperature 0.2:
```bash
python ProteinMPNN/protein_mpnn_run.py --pdb_path_chains B --sampling_temp 0.2 --num_seq_per_target N ...
python ProteinMPNN/protein_mpnn_run.py --use_soluble_model  --sampling_temp 0.2 ...   # SolubleMPNN
```
Output parsed + per-backbone capped → `mpnn_out/filtered.csv` → 2,481 candidate pool.
> **Tip — SolubleMPNN folds better as a monomer** (median pLDDT 93.0 vs 90.3; 90% vs 84%
> clear the gate); it biases to soluble/stable sequences. Useful upstream, but not yet
> evidence of better *binding* — that is decided at the complex stage.

### 02 · localColabFold monomer foldability → 2,162 pass  *(ROOT env, GPU)*
```bash
colabfold_batch --num-recycle 1 --num-models 1 <seqs_fasta_dir> <out_dir>   # pLDDT >= 70 gate
```
> **Tip — monomer foldability barely discriminates.** pLDDT/pTM are correlated (r 0.88–0.91)
> and MPNN sequences fold well as monomers almost by construction (median 91.7; 93.5% ≥ 70).
> Keep the gate **loose** (pLDDT ≥ 70) — over-filtering here discards good binders for no signal.

### 03 · Redundancy + diversity down-sample → 1,050  *(USER venv)*
Dedup at 90% identity (removed **0** here), then diversity down-sample at **70% identity**,
keeping the best-**folding** representative per cluster (needs the pLDDT/pTM from stage 02).
> **Tip — cheap filters first, but check they bite.** The textbook order runs
> dedup/liability/novelty before any GPU fold; confirm they actually remove something.
> 90% dedup removed nothing here (diverse backbones, low MPNN temp), so re-ordering would
> have saved zero folds. The 70% diversity step is deliberately kept *after* folding.

### 04 · localColabFold complex co-fold → 1,050 complexes  *(ROOT env, GPU)*
```bash
bash scripts/run_complex.sh    # recycle 1, 1 model, reused target MSA, 4-chunk resumable loop
```
See the ColabFold time-savers and batching tips above — this is the stage they matter most.

### 05 · Interface scoring  *(USER venv, numpy 1.26.4, CPU)* — **runnable**
```bash
python scripts/ipsae.py       <pae.json> <complex.pdb> 10 10   # -> ipSAE_min
python scripts/score_iface.py <cand>     <complex.pdb>         # -> sc_DockQ, BSA, contacts
python scripts/bench_ddg.py   ...                              # PyRosetta ddG/CMS (constrained full relax)
```
> **Tip — ipSAE is min-of-both-directions.** ipSAE_min = the minimum of the A→B and B→A
> `asym` rows for **one** structure (not a min across models). Pure-PAE, so invariant to
> relaxation — safe to compute once on the raw fold.
>
> **Tip — relax the backbone before ΔΔG.** Side-chain-only relax on recycle-1 structures
> invents clash blow-ups (ΔΔG in the thousands). Use a coordinate-constrained `FastRelax`
> with the backbone free. ipSAE/sc_DockQ barely move; the energy terms become usable.

### 06 · Cloud multi-predictor consensus → 43 pose-PASS  *(USER venv + keys)*
Fold each shortlisted design with three independent cloud co-folders; keep only designs where
**all three** reproduce the pose at `sc_DockQ ≥ 0.23`.
**All three clients are implemented and tested** under `scripts/cofold/`. Each writes
`cand,ipsae_min,sc_dockq` (ipSAE from PAE via `ipsae.py`, sc_DockQ from the predicted
structure via `score_iface.py`) and is resumable.

**Boltz-2** (`boltz2_run.py`) — needs a **live** key; runs under numpy 1.26.4:
```bash
python scripts/cofold/boltz2_run.py --in results/final_selection.csv \
    --out $SCORING/boltz_ipsae_scdockq.csv --scoring $SCORING --scripts scripts --key $BOLTZ_API_KEY
```
`POST api.boltz.bio/.../structure-and-binding`, header `x-api-key`, model `boltz-2.1`; the
archive has `sample_0_predicted_structure.cif` + `sample_0_pae.npz`. Boltz diffusion sampling
is **stochastic** (`num_samples=1`, random seed) — values move a little run-to-run (tested
cand00073: ipSAE 0.90 / sc_DockQ 0.91 vs stored 0.92 / 0.89).

**ESMFold2-Full** (`esmfold2_full.py`) — **two passes**, because the esm SDK needs numpy 2.x
and DockQ needs numpy <2 (they can't coexist):
```bash
pip install -q numpy==2.5.2          # FOLD pass (esm SDK)
python scripts/cofold/esmfold2_full.py fold  --in results/final_selection.csv \
    --outdir /tmp/esm_full_struct --results $SCORING/esmfold2_full_results.csv --token $BIOHUB_TOKEN
pip install -q numpy==1.26.4         # SCORE pass (DockQ)
python scripts/cofold/esmfold2_full.py score --outdir /tmp/esm_full_struct \
    --out $SCORING/esm_full_ipsae_scdockq.csv --scoring $SCORING --scripts scripts
```
Biohub Forge, free + PAE-capable. Tested cand00725: reproduces stored ipSAE 0.8309 / sc_DockQ 0.912 exactly.

**Protenix v2 / JapanFold** (`protenix_run.py`) — free, no key, **rate-limited** so it batches:
```bash
python scripts/cofold/protenix_run.py --in results/final_selection.csv \
    --out $SCORING/protenix_results.csv --scoring $SCORING --scripts scripts \
    --batch 6 --space 3.5 --cooldown 150
```
Submits a small batch (`--batch`, 5–8), spaces submits `--space` s apart, waits for the batch,
then cools down `--cooldown` (~2–3 min) before the next. No PAE (iptm + sc_DockQ only). Tested
cand00073/cand00725: reproduces stored rows exactly.
>
> **Tip — three opinions beat two.** Boltz-2 + ESMFold2-Full agreed on 46 poses; adding
> Protenix rejected 3 more (high Boltz/ESM sc_DockQ but Protenix < 0.23). A cheap third
> co-folder is the highest-value filter in the pipeline.
>
> **Tip — know your free co-folders.** ESMFold2-Full (Biohub) is free + PAE-capable
> (ipSAE works). Protenix (JapanFold) is free but exposes no PAE (iptm + sc_DockQ only).
> Boltz-2 **test** keys return a MOCK fixture — only **live** keys fold real sequences.

### 07 · Rank designs  *(USER venv, numpy 1.26.4)* — **runnable**
```bash
python scripts/analysis/rank_designs.py --scoring $SCORING              # writes final_ranking.csv
python scripts/analysis/rank_designs.py --scoring $SCORING --verify     # check vs existing
```
Composite z = **4·z(ipSAE) + 1·z(sc_DockQ)** (the dataset's 4:1 weighting) over the
pose-passing designs; sort → rank → top-30.

### 08 · Build the interactive report  *(USER venv)*
```bash
python scripts/analysis/novelty_and_extract.py   # novelty vs dataset + extract complexes
python scripts/analysis/slim_report_data.py      # -> report_slim.json (top-30 structures)
python scripts/analysis/align_top30.py           # superpose all top-30 onto the PD-L1 frame
python scripts/analysis/build_report.py          # inject into report.template.html -> report.html
```
Prebuilt report is committed at [`report/index.html`](report/index.html) (open in a browser;
3Dmol.js loads from CDN). These analysis scripts currently use absolute scratch paths — edit
the path constants at the top before re-running elsewhere.

---

## A note for anyone driving this from Windows / PowerShell
`<`, `|`, `$()` in a `wsl … --` command string are parsed by PowerShell first and break the
call. Write the logic to a `.sh`/`.py` file and run it via its `/mnt/c/…` path — far more
reliable than escaping. (All the scripts here are designed to be run that way or from inside WSL.)
