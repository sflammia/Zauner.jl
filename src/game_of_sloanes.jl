export short_fat_matrix, gos_vector, export_gos_file

@doc """
    short_fat_matrix(v::Vector{Complex{BigFloat}})

Takes fiducial vector (length d) as input.
Returns the associated short, fat matrix (d x d^2) given by the Weyl-Heisenberg orbit of v.
"""
function short_fat_matrix(v::Vector{Complex{BigFloat}})
    d = length(v)
    sfm = typeof(v)[]
    for p1 = 0:(d-1)
        for p2 = 0:(d-1)
            push!(sfm, wh(p1,p2,d)*v)
        end
    end
    sfm
end

@doc """
    gos_vector(v::Vector{Complex{BigFloat}})

Takes fiducial vector as input.
Returns flat vector of real parts of vectors in the short, fat matrix followed by the imaginary parts of the vectors in the short, fat matrix.
"""
function gos_vector(v::Vector{Complex{BigFloat}})
    d = length(v)
    sfm = short_fat_matrix(v)
    flatsfm = collect(Iterators.flatten(sfm))
    [real.(flatsfm) ; imag.(flatsfm)]
end

@doc """
    export_gos_file(v::Vector{Complex{BigFloat}}, filename::String)

Takes fiducial vector and a filename as input.
Writes to the filename an encoding of the SIC compatible with the Game of Sloanes database of frames.
Returns the number of characters written.
"""
function export_gos_file(v::Vector{Complex{BigFloat}}, filename::String)
    gos = string(replace(string(gos_vector(v))[10:end-1], ", " => "\n"), "\n")
    write(filename, gos)
end
