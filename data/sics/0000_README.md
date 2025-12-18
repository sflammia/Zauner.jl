This folder contains SIC fiducial vectors computed to floating point precision using `Zauner.jl`.

Each file follows the format
    `k_d_a_b_c.f64`
where
 - `k` is the padded four-digit number in the arbitrary ordering defined by the `dq` function;
 - `d` is the dimension of the SIC;
 - `a`, `b`, `c` are the defining parameters of the quadratic form $Q(x) = ax^2 + bx + c$;
 In this way, $(d,Q)$ is the admissible pair for the associated ghost which generated the SIC via the `necromancy` algorithm.

The association between the SIC fiducial vector and the admissible tuple $(d,Q)$ which labels the original ghost is _noncanonical_ (when the class number is > 1). This is because we have to choose an explicit Galois automorphism to define the mapping from the ghost to the SIC, and this choice is arbitrary (based on our current understanding). This choice could even change in future versions of `Zauner.jl`, so please be careful. The only canonical association possible at this time is across the collection of all SICs in the same multiplet (i.e., pairs of admissible tuples $(d,Q_1)$ and $(d,Q_2)$ where $Q_1$ and $Q_2$ have the same discriminant and conductor, but are in distinct classes).
