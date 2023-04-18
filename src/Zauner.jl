module Zauner

using Reexport, GaussQuadrature, QuadGK, Documenter#, LaTeXStrings
@reexport using Hecke
import Hecke.conductor, Hecke.reduction, Hecke.cycle, Hecke.stabilizer
export reduction, conductor, cycle, stabilizer

include("algebraic.jl")
include("analytic.jl")
include("double_sine.jl")
include("quadform.jl")
include("sl2z.jl")
include("utils.jl")

end # module