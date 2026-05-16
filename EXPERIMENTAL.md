# Experimental branch

The `experiments/ga-exploration` branch holds an in-progress rework of the
Ne-estimation pipeline. It is **not** part of the upstream GONE2 release —
results from it should be considered exploratory, and may diverge from the
master branch.

To try it:

```bash
git checkout experiments/ga-exploration
make gone
# or, for the live ncurses interface:
make gone-ncurses
```

Switch back to the upstream-compatible code with `git checkout master`.

## Main differences vs `master`

- **Genetic algorithm**: the experimental branch runs three GA configurations
  in sequence (truncated quadratic + heavy kicks, L2, L1 + heavy kicks) inside
  a single binary, keeps the pass with the lowest residual, and uses a
  two-stage kick schedule that explores broadly for the first ~40 % of
  generations and then refines.
- **`gone-ncurses` interface**: live chart of the Ne curve as the GA runs,
  with a top-right score box showing the current best residual. The
  experimental build wires this to the combo runtime above.

## Visual comparison

Side-by-side runs of master and the experimental branch on the simulation
suite — including bottlenecks, drops, growth, and constant-Ne cases — are
written up in [docs/branch-comparison.md](docs/branch-comparison.md), with
the rendered charts under [img/](img/).
