```@setup Z
using Zauner
```

# Zauner

`Zauner` is a Julia module for computing properties of SICs (or SIC-POVMs) and related objects such as rank-``r`` versions of SICs (``r``-SICs), objects that are galois conjugates of these called *ghosts*, and certain Stark units over real quadratic fields.
The properties of ghosts, SICs, and Stark units computed by `Zauner` are, in most cases, only known to hold conjecturally.
This package is a systematic way to test Zauner's conjecture and certain Stark conjectures, generate new conjectures, and sharpen (or refute!) existing ones.

!!! warning
    This package is still in its infancy.
    Expect breaking changes in upcoming versions.

## Getting Started

Here's how to get started.
```julia
using Zauner
```

With Julia 1.10 or later, you should be prompted to install `Zauner` if you haven't done so already.
If you're running an older version of Julia, you can first install `Zauner` by running `import Pkg; Pkg.add("Zauner")`.

## Admissible Tuples

The ghosts constructed by `Zauner` are parameterized by an `AdmissibleTuple`.
This can be specified in one of several ways.
A rank-1 ghost in dimension `d > 3` is specified up to equivalence by `d` and a suitable quadratic form `Q`.
Thus, we have
```@repl Z
d = 7
Q = binary_quadratic_form(2,-4,1)
t = AdmissibleTuple(d,Q)
```
We see that the `AdmissibleTuple` type precomputes several associated data from `d,Q`.
```@repl Z
t.K # associated field
t.h # class number / size of multiplet
(t.j, t.m) # dimension grid positions
t.D # fundamental discriminant
```
The complete set of fields is given in the documentation for `AdmissibleTuple`.

If no quadratic form is given, then `AdmissibleTuple` defaults to the principle form.
```@repl Z
t = AdmissibleTuple(7)
t.Q # form
```

An `AdmissibleTuple` can be used to do explicit computations of ghosts.
From the previous tuple, we have the following ghost (reduced from default `BigFloat` precision):
```@repl Z
v = ComplexF64.(ghost(t))
```
Note that ghosts are normalized so that the first coordinate is 1.
To get the properly normalized ghost overlaps and validate the solution, one can do the following to compute ghost overlaps:
```@repl Z
u = circshift(reverse(v),1); # parity-reversed ghost
nu_plus  = real.( [u' * wh( [ p; q], v) / (u'v) for p=0:d-1, q=0:d-1] )
nu_minus = real.( [u' * wh( [-p;-q], v) / (u'v) for p=0:d-1, q=0:d-1] )
(d + 1) .* nu_plus .* nu_minus
```
Here `wh(p,v)` acts a Weyl-Heisenberg displacement operator indexed by `p` onto `v`, i.e., in quantum notation it returns $D_{\boldsymbol{p}}|v\rangle$.

It can also be used with `necromancy` to attempt to reconstruct a 1-SIC, which we again round down from `BigFloat` precision.
```@repl Z
ψ = ComplexF64.(necromancy(t))
```

`Zauner` provides numerical validation by measuring how far the returned objects are from satisfing the overlap conditions or minimizing the "pointwise" frame potential conditions discussed below.
```@repl Z
ghost_overlap_test(v)
sic_overlap_test(ψ)
ghost_frame_test(v)
sic_frame_test(ψ)
```

!!! warning
    Rank ``r>1`` has very limited support at the moment and is essentially confined to the definition of an `AdmissibleTuple`.
    Many other functions have limited domains at present.
    For example, we do not presently support necromancy on ``F_a`` orbits.
    Please file an issue on GitHub if there is something you wish to see supported in future updates.

## Citing `Zauner`

How to cite us:
(**add citation info for the GitHub repo and for the paper**)
```
@misc{AFK2024,
	archiveprefix = {arXiv},
	author = {Appleby, Marcus and Flammia, Steven T. and Kopp, Gene S.},
	eprint = {2407.00001},
	howpublished = {arXiv:2407.00001},
	primaryclass = {math.NT},
	title = {A constructive approach to Zauner's conjecture via the Stark conjectures},
	year = {2024}
}
```
