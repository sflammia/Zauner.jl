# This script generates the data tables in the appendix in latex format

parent_dir_name = dirname(@__DIR__)
table_name = "data.tex"
file_name = parent_dir_name * "/src/" * table_name


# pretty print the cycle decomposition of v
function _print_cycle_decomp(io, v::AbstractVector)
    if v == []
        return print(io, " C_{1} ")
    end
    unique_elements = Dict{Int,Int}()
    for item in v
        if haskey(unique_elements, item)
            unique_elements[item] += 1
        else
            unique_elements[item] = 1
        end
    end
    c_strings = []
    for key in sort(collect(keys(unique_elements)))
        value = unique_elements[key]
        if value > 1
            push!(c_strings, "C_{$key}^{$value}")
        else
            push!(c_strings, "C_{$key}")
        end
    end
    print(io, join(c_strings, "\\times{}"))
    print(io, " ")
end


open(file_name, "w") do io
    println(io, "\\begin{center}")
    println(io, "\\begin{longtable}{tttttttttt}")
    println(io, "\\toprule")
    println(io, "d & \\Delta_0 & f & h & \\Cl(\\mathcal{O}_f) & \\Gal(E/H) & Q & L^n & \\text{a.u.} & \\ell \\\\")

    d = 1
    q = 0
    for j = 1:614
        F = AdmissibleTuple(dq(j))
        au, L = is_antiunitary_with_generator(F)
        if !(au)
            L = F.L
        end
        print(io, (F.d == d ? "" : "\\midrule\n"))
        print(io, (F.d == d ? "" : F.d), "&") # dimension d
        print(io, (F.d == d ? "" : F.D), "&") # fundamental discriminant Δ
        print(io, (F.d == d && F.q == q ? "" : F.q), "&") # conductor f
        if !(F.d == d && F.q == q)
            cgp = class_group_structure(F.q^2 * F.D)
            print(io, " $(prod(cgp)) & ") # class number
            _print_cycle_decomp(io, cgp)
        else
            print(io, " &")
        end
        print(io, " & ")
        if !(F.d == d && F.q == q)
            ords = galois_group_structure(F)
            _print_cycle_decomp(io, ords) # galois group
        end
        print(io, "&\\langle", F.Q.a, ",", F.Q.b, ",", F.Q.c, "\\rangle&")
        print(io, "\\smt{", L[1, 1], "&", L[1, 2], "\\\\", L[2, 1], "&", L[2, 2])
        print(io, "}^{", (au ? 2 : 1) * F.k, "}&")
        print(io, (au ? "\\text{Y}&" : "&")) # is anti-unitary?
        println(io, length(psl2word(F.A)), "\\\\")
        d = F.d
        q = F.q
    end

    println(io, "\\bottomrule")
    println(io, "\\end{longtable}")
    println(io, "\\end{center}")
end
