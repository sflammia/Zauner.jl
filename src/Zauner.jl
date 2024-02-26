module Zauner

using LazilyInitializedFields
using Reexport, Documenter#, LaTeXStrings
@reexport using Hecke
using SpecialMatrices: Vandermonde
using AMRVW: roots
using ForwardDiff: jacobian
using GenericFFT: ifft, fft
using QuadGK: quadgk
import Hecke.conductor, Hecke.cycle, Hecke.reduction, Hecke.stabilizer
export conductor, cycle, reduction, stabilizer

include("utils.jl")
include("algebraic.jl")
include("analytic.jl")
include("double_sine.jl")
include("quadform.jl")
include("sl2z.jl")
include("precision_bump.jl")
include("galois.jl")

end # module