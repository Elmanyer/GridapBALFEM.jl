# SNELLIUS_ROME_LAUNCH_CONFIGS.md — how the `rome` partition charges, and how to size a job

> **Tracked deliberately.** This governs what every cluster job *costs*, so it lives beside the
> launchers in `run/` rather than in the gitignored `markdown_files/`.
>
> Cluster: **SURF Snellius**, partition **`rome`**. Enforced by Slurm **cgroups**: what you request
> is exactly what you get, and exactly what you are billed for.

---

## 1. The one rule that decides your bill

A `rome` node is **128 cores / 224 GiB**, and it is divided into **eighths**. You cannot rent less
than an eighth, and you cannot rent a fraction that is not a multiple of one.

```
1/8 node  =  16 cores  and  28 GiB
```

> ### ⚠ THE BILLING RULE
> **Your reservation share is `max(cores/128, memory/224)`, rounded UP to the next eighth.**
> Cores and memory are evaluated *separately* and **the larger fraction wins**. If your memory
> request exceeds the proportional share of your core request, the allocation **steps up to the next
> valid boundary** and you pay for cores you never asked for.

The consequence worth internalising: **1/8 of a node gives you 28 GiB for 16 cores, i.e. 1.75 GiB
per core.** That ratio — not the core count — is the thing to size against.

```
mem-per-cpu  <  1.75 GiB   →  CORES decide the bill   (memory is free capacity you are wasting)
mem-per-cpu  =  1.75 GiB   →  perfectly balanced      (both fractions equal)
mem-per-cpu  >  1.75 GiB   →  MEMORY decides the bill (cores are idle and you pay for them)
```

**BALFE-M needs ≈4 GiB per rank** (measured: baseline ≈1.4 GB before solving anything, peak
3633 MB/rank on the small-domain suite; a job at the 2 GB/core node default was OOM-killed). So this
solver is **permanently in the third regime — memory always decides the bill.** Every sizing
decision below follows from that single fact.

---

## 2. The valid allocations

| fraction | cores | memory | ranks @ **4 GiB** (BALFE-M) | cores left idle |
|---|---|---|---|---|
| 1/8 | 16 | 28 GiB | **7** | 9 |
| 2/8 (1/4) | 32 | 56 GiB | **14** | 18 |
| 3/8 | 48 | 84 GiB | **21** | 27 |
| 4/8 (1/2) | 64 | 112 GiB | **28** | 36 |
| 5/8 | 80 | 140 GiB | **35** | 45 |
| 6/8 (3/4) | 96 | 168 GiB | **42** | 54 |
| 7/8 | 112 | 196 GiB | **49** | 63 |
| 8/8 (full) | 128 | 224 GiB | **56** | 72 |

**`ranks = 7 × k` for tier `k/8`.** At 4 GiB/rank you can never use more than 43.75 % of the cores in
your own allocation — the idle-core column is not waste you can recover, it is the price of the
memory footprint. Do not try to "fill" it by adding ranks; that just steps you up a tier.

### Worked examples — the same job, sized three ways

| request | memory | core frac | mem frac | **billed** | verdict |
|---|---|---|---|---|---|
| 32 ranks @ 4G | 128 GiB | 2/8 | **4.57/8** | **5/8** | ⚠ pays for 80 cores, uses 32 |
| 28 ranks @ 4G | 112 GiB | 1.75/8 | **4.00/8** | **4/8** | ✅ exact fit, 20 % cheaper |
| 48 ranks @ 4G | 192 GiB | 3/8 | **6.86/8** | **7/8** | ⚠ 6 GiB from a full node |
| 42 ranks @ 4G | 168 GiB | 2.63/8 | **6.00/8** | **6/8** | ✅ exact fit |
| 40 ranks @ 5G | 200 GiB | 2.5/8 | **7.14/8** | **8/8** | ⛔ a whole node for 40 ranks |

⚠ **`--mem-per-cpu=5G` is the worst of both worlds** and is not used anywhere in this repository:
it buys nothing the solver needs and pushes almost every rank count over a tier boundary.
**Use 4G.**

---

## 3. Sizing procedure

1. **Pick the rank count the physics wants** — enough DOFs per rank to be worth decomposing
   (`CLAUDE.md` rule 23), cells near-isotropic (rule 22).
2. **Compute the tier:** `k = ceil(4 × ranks / 28)`.
3. **Snap to the tier's exact rank count, `7k`.** If the physics number falls *between* two tiers,
   **take the LOWER tier** — pay less and run slightly fewer ranks, rather than pay a whole extra
   eighth for a handful of ranks.
4. **Factor it into `PX × PY`** (`7k` factors as `7 × k`, which is why `7×4 = 28` and `7×6 = 42`
   drop straight in where `8×4 = 32` and `8×6 = 48` used to sit).
5. **Multi-node:** the rule applies **per node**. Request `ntasks-per-node = 7k` and pay `k/8` on
   each node.

```bash
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G        # 112 GiB = exactly 4/8 of a rome node
export BALFEM_PX=7
export BALFEM_PY=4              # 7*4 = 28
```

---

## 4. Traps

* **A rome node advertises 256 GB but Slurm allocates ~224 GiB.** `64 × 4G = 256 GiB` is refused at
  submit time. The table above uses the real 224.
* **`--mem` and `--mem-per-cpu` are exclusive.** This repo uses `--mem-per-cpu` throughout so the
  arithmetic stays per-rank.
* **One node beats two, when the problem fits.** All six production launchers were moved from
  `2 × 42 ranks` (6/8 on each of two nodes) to `1 × 42` (6/8 on one) on 2026-09-07.
  ⚠ **Be precise about why this is cheaper.** It halves the *reservation*, but billing is
  fraction × wall-time, and halving the ranks roughly doubles the wall-time — so on a perfectly
  scaling problem the node-hours would come out **the same**. It wins because scaling is *not*
  perfect: at 84 ranks each subdomain is small, the GMRES communication share is larger, and
  parallel efficiency is lower than at 42. Fewer, better-loaded ranks therefore cost fewer
  node-hours for the same work — and a single-node job also queues sooner and avoids inter-node
  communication entirely. **Re-measure if a case grows; this is an efficiency argument, not an
  identity.**
* **Fewer ranks is often cheaper *and* faster.** Beyond the point where each rank holds a meaningful
  share, adding ranks raises the communication cost *and* the bill. Spend spare capacity on more
  *cases*, not on decomposing one case further (`CLAUDE.md` rule 23).
* **The sysimage build is a single task with many threads**, so it is sized by
  `cpus-per-task × mem-per-cpu`, not by rank count — the same billing rule applies.
* **`run_blue.sh` targets DelftBlue, not Snellius.** Its `3900M` is correct for that machine and
  none of this applies to it.
