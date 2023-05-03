module Zauner

using Reexport, Documenter#, LaTeXStrings
@reexport using Hecke
import QuadGK: quadgk
import GenericFFT: ifft, fft
import Hecke.conductor, Hecke.reduction, Hecke.cycle, Hecke.stabilizer
export reduction, conductor, cycle, stabilizer

include("algebraic.jl")
include("analytic.jl")
include("double_sine.jl")
include("quadform.jl")
include("sl2z.jl")
include("utils.jl")

end # module