module Zauner

using Reexport, Documenter
@reexport using Hecke
using LazilyInitializedFields
using LinearAlgebra
using SpecialMatrices: Vandermonde
using AMRVW: roots
using ForwardDiff: jacobian
using GenericFFT: ifft, fft
using QuadGK: quadgk
using Convex
using SCS
import Hecke.conductor, Hecke.cycle, Hecke.reduction, Hecke.stabilizer
export conductor, cycle, reduction, stabilizer

include("utils.jl")
include("algebraic.jl")
include("hj.jl")
include("analytic.jl")
include("double_sine.jl")
include("quadform.jl")
include("sl2z.jl")
include("precision_bump.jl")
include("integer_relations.jl")
include("galois.jl")
include("matrix_completion.jl")
include("invariants.jl")

end # module
