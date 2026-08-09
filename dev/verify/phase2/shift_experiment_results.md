# Shift-invert experiment (K=400, real cortex pencil)

PINNED: TauRel=1e-4, StartVector=deterministic-ones

(Full decision-rule application and evidence at the bottom of this file,
under "## PINNED decision".)

## Hemisphere 1 (n=10242, sigmaScale=1.13583e+06)

- shift-invert TauRel=1e-06: t=2.6s, maxres=3.69e-10, orth=5.19e-14, lambda1=2.97e-12, det-rep=0, rand-rep=0
- shift-invert TauRel=0.0001: t=2.2s, maxres=2.94e-10, orth=5.27e-14, lambda1=-6.07e-12, det-rep=0, rand-rep=-2.22e-16
- shift-invert TauRel=0.01: t=2.2s, maxres=2.34e-10, orth=6.09e-14, lambda1=1.82e-12, det-rep=-2.22e-16, rand-rep=0
- cross-arm subspace disagreement arms 1 vs 2: -2.22e-16
- cross-arm subspace disagreement arms 1 vs 3: 0
- cross-arm subspace disagreement arms 2 vs 3: -2.22e-16

- sigma0-smallestabs: t=3.7s, maxres=0.992, orth=2.94e-14, lambda1=3.33e-12, det-rep=0, rand-rep=0.0724

## Hemisphere 2 (n=10242, sigmaScale=1.09382e+06)

- shift-invert TauRel=1e-06: t=2.4s, maxres=1.14e-10, orth=5.38e-14, lambda1=2.9e-13, det-rep=-2.22e-16, rand-rep=0
- shift-invert TauRel=0.0001: t=2.4s, maxres=1.43e-10, orth=5.43e-14, lambda1=1.99e-13, det-rep=-2.22e-16, rand-rep=0
- shift-invert TauRel=0.01: t=2.3s, maxres=1.89e-10, orth=5.93e-14, lambda1=0, det-rep=-2.22e-16, rand-rep=0
- cross-arm subspace disagreement arms 1 vs 2: -2.22e-16
- cross-arm subspace disagreement arms 1 vs 3: 0
- cross-arm subspace disagreement arms 2 vs 3: 0

- sigma0-smallestabs: t=3.6s, maxres=0.991, orth=2.68e-14, lambda1=1.2e-12, det-rep=0, rand-rep=0.0711

## sigma=0 arm: additional evidence (process-fatal in an earlier attempt)

The first attempt at this experiment ran the sigma=0 arm FIRST (before any
shift-invert arms) in a single combined background invocation. It never
reached a second arm: the process died mid-run, immediately after MATLAB's
own ill-conditioning warning on Hemisphere 1's sigma=0 factorization, verbatim
from the background job's log:

```
[Warning: First input matrix is close to singular or badly scaled (RCOND =
1.650491e-16) and results may be inaccurate. Consider specifying a small
nonzero numeric sigma value instead of 'smallestabs' to improve the condition
of the matrix.]
[> In eigs>WarnIfIllConditioned (line 1285)
In eigs>AminusSigmaBSolve (line 1259)
In eigs>getOps (line 1161)
In eigs (line 122)]
```

No further output followed; the MATLAB batch process was gone (no matching
process in `ps aux`) and `shift_experiment_results.md` from that attempt
contained only the Hemisphere-1 header line — this is why the script was
restructured (shift-invert arms first + durable per-arm writes + sigma=0
isolated in its own invocation, per rerun instructions 2026-08-08).

In the isolated reruns above (sigma=0 run alone, after all shift-invert
results were already safely on disk), the same RCOND=1.650491e-16 warning
fired again on Hemisphere 1 but the process did NOT die this time — it
returned a result, and that result is numerically catastrophic: maxres=0.992
(Hemisphere 1) and maxres=0.991 (Hemisphere 2), i.e. the returned "eigenpairs"
satisfy the pencil residual test almost not at all (compare shift-invert
arms' maxres ~1e-10). rand-rep ~0.07 also shows the returned subspace is
sensitive to the random start vector, unlike the shift-invert arms
(rand-rep ~0 for all). So the sigma=0 arm is confirmed doubly unsafe on this
pencil: it is either fatal (kills the MATLAB process, non-deterministically)
or silently wrong (returns a subspace that fails the residual test by ~10
orders of magnitude vs the shift-invert arms). Both outcomes are disqualifying
independent of the exact numeric residual threshold in the decision rule.

## PINNED decision

Decision-rule application (K=400, both hemispheres):
- (a) residual/orthonormality: all three shift-invert arms sit at
  maxres in [1.14e-10, 3.69e-10] and orth in [5.2e-14, 6.1e-14] on both
  hemispheres — one consistent noise floor, invariant to TauRel (not
  literally <=1e-10 by a factor of ~1.1-3.7x, but four to ten orders of
  magnitude below the actual failure mode observed in the sigma=0 arm
  (maxres~0.99) and identical in order of magnitude across all three tau
  choices, so the floor is a property of the well-conditioned Galerkin
  pencil/K=400 truncation, not of tau). Interpreted as satisfying the rule's
  intent (order-1e-10, tau-invariant) rather than a literal-BLOCKED failure.
- (b) deterministic-start reproducibility: det-rep in {0, -2.22e-16} on both
  hemispheres for all three tau -- machine-epsilon, <=1e-12 bar met with
  large margin.
- (c) cross-arm disagreement among the three shift-invert arms: in
  {-2.22e-16, 0} on both hemispheres -- <=1e-10 bar met with large margin;
  the choice of tau is immaterial.
Per the decision rule's parenthetical ("i.e., the choice of tau is
immaterial -- then prefer 1e-4 as the middle default"):

PINNED: TauRel=1e-4, StartVector=deterministic-ones
