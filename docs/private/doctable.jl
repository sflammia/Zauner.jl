# This script generates the table for the documentation
# using Zauner

parent_dir_name = dirname(@__DIR__)
table_name = "tables.md"
file_name = parent_dir_name * "/src/" * table_name

# pretty print the cycle decomposition of v
function _print_cycle_decomp(io, v::AbstractVector)
    if v == []
        return print(io, "\$ C_{1} \$")
    end
    print(io, "\$")
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
    print(io, "\$ ")
end



preamble = "# Admissible tuple data tables\n
This table contains algebraic data for each inequivalent admissible tuple in dimensions 4--35, comprising 100 total tuples.
This list is conjecturally complete for all Weyl-Heisenberg covariant SICs in these dimensions.
Each admissible tuple is specified by a dimension \$d\$ and an integer binary quadratic form \$Q\$ as follows.
First factorize \$(d+1)(d-3)=s^2\\Delta_0\$ where \$\\Delta_0\$ is a fundamental discriminant.
Then \$ (d,Q) \$ gives an admissible tuple if \$\\mathrm{disc}(Q) = f^2\\Delta_0\$ where \$f\$ divides \$s\$.
The other columns can be computed from these data, but they may be difficult to compute,
for example requiring integer factoring or finding a fundamental unit.
The column \$\\Delta_0\$ contains the fundamental discriminant of \$Q\$ and \$h\$ is the order of the class group \$\\mathrm{Cl}(\\mathcal{O}_f)\$, given in the next two columns respectively.
The Galois group \$\\mathrm{Gal}(E/K)\$ of the field containing the overlaps over \$K=\\mathbb{Q}(\\sqrt{\\Delta_0})\$ is given in the next column.
As both the class group and the Galois group are finite and abelian, we give the canonical decomposition into cyclic groups \$C_k\$ of order \$k\$.
For the special case that the tuple has so-called \$F_a\$ symmetry, we have not yet worked out the Galois groups, so we mark these entries as tbd.
The column \$L^n\$ contains a generator \$L\$ of the stability group of \$Q\$ in \$\\mathrm{GL}_2(\\mathbb{Z})\$ and its order \$n\$ in \$\\mathrm{GL}_2(\\mathbb{Z}/\\bar{d})\$;
that is, treating \$Q\$ as a symmetric matrix we have \$L^T Q L = \\det(L) Q\$ and \$L^n = 1\\ (\\bmod\\ \\bar{d})\$.
If the tuple has antiunitary symmetry, we denote this with a Y in the a.u. column.
Finally, \$\\ell\$ is the length of the word expansion of \$L^n\$ using the Hirzebruch-Jung (negative regular) reduction into the standard
(\$S\$ and \$T\$) generators of \$\\mathrm{SL}_2(\\mathbb{Z})\$.
This is one measure of the complexity of constructing the actual fiducial vector for that input.
The \$Q\$ in this list were selected among class representatives to minimize this complexity, although this choice is not generally unique.
The data here are sufficient to compute a ghost fiducial in each class,
but to fully specify a SIC,
one must additionally choose a sign-switching Galois automorphism \$\\sqrt{\\Delta_0}\\to-\\sqrt{\\Delta_0}\$ over an appropriate field extension of \$K\$.
"

open(file_name, "w") do io
    println(io, preamble)
    println(io, "| \$d\$ | \$\\Delta_0\$ | \$f\$ | \$h\$ | \$\\mathrm{Cl}(\\mathcal{O}_f)\$ | \$\\mathrm{Gal}(E/K)\$ | \$Q\$ | \$L^n\$ | \$\\text{a.u.}\$ | \$\\ell\$ |")
    println(io, "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    d = 1
    q = 0
    # for j=1:282 # 4 ≤ d ≤ 64
    for j = 1:100 # 4 ≤ d ≤ 35
        # for j=1:10 # testing
        F = AdmissibleTuple(dq(j))
        au, L = is_antiunitary_with_generator(F)
        if !(au)
            L = F.L
        end
        print(io, "|")
        print(io, (F.d == d ? " |" : " \$$(F.d)\$ |")) # dimension d
        print(io, (F.d == d ? " |" : " \$$(F.D)\$ |")) # fundamental discriminant Δ
        print(io, (F.d == d && F.q == q ? " |" : " \$$(F.q)\$ |")) # conductor f
        if !(F.d == d && F.q == q)
            cgp = class_group_structure(F.q^2 * F.D)
            print(io, " \$$(prod(cgp))\$ | ") # class number
            _print_cycle_decomp(io, cgp)
        else
            print(io, " |")
        end
        print(io, " | ")
        if !(F.d == d && F.q == q)
            if has_fa_symmetry(F)
                print(io, " tbd ") # galois group
            else
                ords = galois_group_structure(F)
                _print_cycle_decomp(io, ords) # galois group
            end
        end
        print(io, " | \$\\langle", F.Q.a, ",", F.Q.b, ",", F.Q.c, "\\rangle\$|\$")
        print(io, "\\left(\\begin{smallmatrix}", L[1, 1], "&", L[1, 2], "\\\\", L[2, 1], "&", L[2, 2])
        print(io, "\\end{smallmatrix}\\right)^{", 2^au * F.k, "}\$|")
        print(io, (au ? "\$\\text{Y}\$" : ""), "|\$") # is anti-unitary?
        println(io, length(psl2word(F.A)), "\$|")
        d = F.d
        q = F.q
    end
end
