# Lessons learned — PD-L1 binder pipeline (for final report)

## 1. Filter ordering: cheap (sequence-level) before expensive (GPU)
**Principle.** Sequence-only filters — redundancy dedup, liability (Cys parity, homopolymer),
novelty (MMseqs2) — cost milliseconds and should run BEFORE any GPU folding
(monomer ~11 s/seq, complex much more). The Anthropic PD-L1 prompt groups all of these under
"Pre-scoring filters — run before any co-folding spend."

**What we actually did.** We ran the monomer-foldability screen on all 2,481 selected
sequences first, then dedup/liability. Ideal order would put dedup+liability+novelty first.

**Nuance (important, keep it honest in the report).**
- Redundancy dedup at the spec threshold (90% identity) removed **0** sequences on this pool
  (52 diverse RFdiffusion backbones + ProteinMPNN temp 0.2 give no near-duplicates at 90%).
  So doing it before the monomer screen would have saved **0** folds — no compute was wasted
  in this run. That was luck, not design.
- The step that genuinely shrinks the pool is **diversity down-sampling at 60-70% identity**
  (0.6 -> 331 clusters, 0.7 -> 1,050). That step should run **AFTER** monomer folding, because
  we want to keep the best-folding representative per cluster, which needs the pLDDT/pTM scores.
  Clustering-then-picking before folding forces a blind pick by MPNN global_score (weaker).

**Corrected optimal pipeline:**
`dedup 90% + liability + novelty (cheap)` -> `monomer-fold` ->
`diversity down-sample using fold scores to pick best-per-cluster` -> `complex-fold` ->
`rank by ipSAE_min + sc_DockQ` -> `novelty vs UniRef90 before final 30`.

## 2. Monomer foldability is a weak discriminator for binders
pLDDT and pTM are strongly correlated (r 0.88-0.91), and MPNN sequences fold well as monomers
almost by construction (median pLDDT 91.7; 93.5% >= 70). The monomer gate is a sanity filter,
not the ranking signal — binder quality is decided at the complex stage (ipSAE_min + sc_DockQ).
Keep the monomer gate loose (pLDDT >= 70, their default), don't over-filter here.

## 3. SolubleMPNN folds better than ProteinMPNN as monomer
Median pLDDT 93.0 vs 90.3; pLDDT>=90 66% vs 52%; monomer-gate pass 90% vs 84%. Expected
(SolubleMPNN biases toward soluble/stable sequences). Not yet evidence of better binding.

## 4. Novelty at small candidate count
UniRef90 DB size is fixed regardless of query count, so "fewer candidates" does not shrink the
install. But with only tens of candidates you can substitute a local UniRef90 host with online
NCBI BLASTp vs nr + local targeted checks (target chains, Ubiquitin P0CG47/48, known-binder
corpus). Their spec requires full UniRef90 only for the FINAL sheet; online substitution is a
disclosed deviation appropriate for a resource-constrained reproduction.
