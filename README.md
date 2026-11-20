# Experiments on support minors with hints

These experiments implement the full-column reduction of Section 4 and the matching-refined Support Minors method of Sections 5.1–5.3. Each trial starts with a planted, rank-\(r\) MinRank instance over a prime field. The mask consists of exactly \(h\) uniformly sampled positions of \(C'\), sampled without replacement. The ten masks for a parameter set are tested on the same public instance, so the experiment illustrates variation due to the hint locations. These toy instances are not estimates of the security of cryptographic parameter sets.

## Summary

All 50 trials met the stated criterion within the tested degrees. The cost in this table is the finite-field operation count for row-reducing the Macaulay matrix at the successful degree only; it does not include preprocessing or unsuccessful lower-degree attempts. The extrema of the degree and cost are computed independently over the ten trials.

| Parameters (q, m, n, K, r) | h | Solved | Smallest b | Smallest cost | Largest b | Largest cost |
|---|---:|---:|---:|---:|---:|---:|
| (7, 4, 5, 5, 2) | 2 | 10/10 | 1 | 413 | 1 | 44,566 |
| (7, 4, 5, 6, 2) | 3 | 10/10 | 1 | 881 | 1 | 35,456 |
| (7, 5, 6, 7, 2) | 3 | 10/10 | 1 | 6,227 | 1 | 307,356 |
| (7, 5, 6, 8, 2) | 4 | 10/10 | 1 | 357 | 1 | 228,564 |
| (7, 5, 6, 7, 3) | 5 | 10/10 | 1 | 2,058 | 2 | 10,425,984 |

The degrees are the *first tested degree for which this implementation recovers and validates a solution*, summarized across the ten trials of each set. Their minimum and maximum need not correspond to the minimum and maximum cost, respectively. A blank result means the stated criterion was not met over the tested degrees; it does not establish that the instance is unsolvable at those or higher degrees. The full trial records and the individual degrees attempted are in the supplied results file.

## Procedure

Write \(C=(I_r\mid C')\) and \(E(x)=\sum_{i=1}^Kx_iM_i\). If \(a\) columns of \(C'\) are completely revealed, impose the \(m\) homogeneous constraints for each such column \(j\):
\[
 E(x)_{k,r+j}-\sum_{i=1}^r C'_{i,j}E(x)_{k,i}=0\quad(1\le k\le m).
\]
The implementation computes their actual rank \(\rho\), parametrizes the nullspace with \(d=K-\rho\) coordinates, substitutes this basis into the public matrices and deletes the fully known columns. Thus it does not assume that the \(am\) equations are independent.

It then subtracts the hinted entries from the remaining nonidentity columns using identity-column operations. For each \(r\)-subset \(T\), it tests whether the normalized support pattern has a perfect matching and keeps its formal minor coordinate if so. It retains a block indexed by an \((r+1)\)-subset \(L\) when some \(L\setminus\{j\}\) survives; each retained block supplies a support-minor equation for each row. The matching is based only on the mask and uses the public revealed values in the normalization; it does not require access to the unknown entries.

For each tested \(b\), the script multiplies the resulting bilinear equations by all \(x\)-monomials of degree \(b-1\) and performs dense finite-field row reduction. With \(s_\Omega\) surviving minor coordinates and \(\ell_\Omega\) retained blocks, the full column count is
\[
N_b=s_\Omega\binom{d+b-1}{b},
\]
and the maximum generated row count is \(m\ell_\Omega\binom{d+b-2}{b-1}\); identically zero equations can make the actual row count smaller. The results report the actual dimensions and rank, rather than assuming the rows are independent.

The operational success criterion is that the Macaulay matrix have a one-dimensional right kernel and that a solution decoded from the kernel validate against the original instance and hints. The identity-minor coordinate is used to extract the \(x\)-ratios. Corank one is sufficient for this linearization procedure, but it is not necessary for solvability. A larger kernel may still contain a valid solution; these experiments do not exhaustively search such kernels. Only degrees \(1\le b<q\) within the configured limit are tested.

## Cost conventions

* **Linear reduction:** `linear_rref_field_ops` counts the additions, subtractions, multiplications, and inversions performed by the finite-field RREF of the full-column constraint matrix. Forming that matrix is timed but not included in this count.
* **Basis substitution and normalization:** `transform_field_ops_model` is a transparent arithmetic model for the implemented computation, counting \(2K\) operations for each transformed matrix entry and two for each normalization update. It includes operations with zero coefficients and is not a count for an optimized algorithm.
* **Matching:** `matching_r_subset_tests` counts the \(r\)-subsets tested, not the permutation checks within each test. The implementation searches permutations and can take up to \(r!\) checks per subset. Unlike the polynomial-time matching procedure described in the draft, this brute-force implementation is intended only for small \(r\). This stage also finds the retained equation blocks; `matching_seconds` measures both tasks.
* **Equation and matrix assembly:** `equations_seconds` gives the time to assemble the support-minor equations. `build_seconds` and `nnz` in the attempts file describe each Macaulay matrix.
* **Macaulay solve:** `elim_field_ops` counts each addition, subtraction, multiplication, and inversion executed in the RREF as one operation, with `pivot_inversions` reported separately. `cumulative_macaulay_elim_field_ops` sums this count over all attempted degrees for a trial. It excludes basis construction, matrix assembly, candidate decoding, and validation.
* **Elapsed times:** the trial file reports measured times for the main phases and the overall run; the attempts file reports matrix construction and row reduction times by degree. The overall time includes decoding and validation, although those tasks do not have a separate timer. Wall-clock measurements depend on hardware, operating-system load, and the Python implementation and may vary between runs.

These are concrete costs for the implementation and hardware used here, rather than the draft’s sparse Wiedemann estimate. In particular, an actual dense elimination count is not the complexity in equation (3) of the paper, and a test count is not a finite-field operation count. Costs depend on the implementation, data ordering and whether all lower degrees have already been attempted. The reported minimum and maximum are empirical over the ten sampled masks, not theoretical bounds over all placements.

## Reproducing the experiment

Use Python 3 with its standard library; no external packages are needed. In the directory containing `experiment_agent.py`, run:

```bash
python experiment_agent.py
```

This generates the 50 seeded trials and writes `experiment_results_agent.json`. Elapsed times will vary on a rerun. To solve a supplied instance, run `python experiment_agent.py example_instance.json`. An input JSON file has keys `q`, `r`, `matrices`, `hints`, and optionally `max_b`. The list `matrices` contains the \(K\) public \(m\times n\) matrices; a hint `[i, j, value]` refers to row `i` and column `j` of \(C'\), both zero-based. The implementation assumes the first \(r\) columns are the identity chart and that a rank-\(r\) solution exists in this chart. The bundled example is a small instance to test the input format.

The raw JSON preserves every record. `trials.csv` has one row per mask, including actual revealed values and preprocessing costs; `macaulay_attempts.csv` has one row per degree tried. A dense matrix exceeding two million entries is not built: that attempt is marked `resource_cap`. For the bundled experiments the maximum is \(b=3\). The command-line default is \(\min(3,q-1)\); an explicit maximum must lie between 1 and \(q-1\).

## Scope and limitations

* The solver targets solutions of **exact rank** \(r\) whose first \(r\) columns are independent. It does not find rank-below-\(r\) solutions or solutions requiring a different identity chart.
* When every column of \(C\) is revealed, the code searches the reduced projective space directly, with a cap of 10,000 candidates; none of the 50 trials takes this branch.
* All fields used here are prime fields; extending the arithmetic to \(\mathbb F_{p^k}\), such as \(\mathbb F_{16}\), requires implementing extension-field operations.
* These instances are planted and conditioned to have an invertible identity chart. They are small examples, not samples from an unspecified production scheme. The rank and success frequencies should not be extrapolated to larger parameters.
* Independent hint placements do not imply independent public instances: each parameter set uses one public instance and ten masks. This design helps isolate the effect of the mask; it does not estimate variation over random instances.
* The Macaulay model relaxes the determinantal relations between the surviving formal minors. Verifying a decoded candidate against the original system and the revealed entries is therefore essential. The rank of the Macaulay matrix and the success criterion are empirical for each run.
* A trial that never reaches the stated criterion is reported as unsuccessful within the tested range, not as an impossible MinRank instance. A different solving algorithm or a larger degree may succeed.

To reproduce a Method 1 baseline on a given input, pass only the hints that lie in fully revealed columns. The solver then performs the full-column reduction and runs the same Support Minors solver on the reduced instance without the partial-column hints. A baseline is not included in the reported table.
