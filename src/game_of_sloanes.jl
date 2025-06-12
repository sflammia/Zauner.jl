export short_fat_matrix, gos_vector

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

function gos_vector(v::Vector{Complex{BigFloat}})
	d = length(v)
	sfm = short_fat_matrix(v)
	flatsfm = collect(Iterators.flatten(sfm))
	[real.(flatsfm) ; imag.(flatsfm)]
end

function export_gos_file(v::Vector{Complex{BigFloat}}, filename::String)
	gos = string(replace(string(gos_vector(v))[10:end-1], ", " => "\n"), "\n")
	write(filename, gos)
end