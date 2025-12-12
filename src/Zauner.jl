module Zauner

using Reexport
ENV["HECKE_PRINT_BANNER"] = false
@reexport using Hecke
using ForwardDiff: jacobian
using GenericFFT: fft, ifft
using LazilyInitializedFields
using LinearAlgebra
using QuadGK: quadgk

import Hecke.conductor, Hecke.cycle, Hecke.reduction, Hecke.stabilizer, Hecke.QuadBin
export conductor, cycle, reduction, stabilizer

include("dQ.jl")
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
include("fourier.jl")
include("invariants.jl")
include("validation.jl")

end # module
