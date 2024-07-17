# This script generates the html version of the data table

parent_dir_name = dirname(@__DIR__)
table_name = "table.html"
file_name = parent_dir_name * "/src/" * table_name

# pretty print the cycle decomposition of v
function _print_cycle_decomp(io, v::AbstractVector)
    if v == []
        return print(io, "<msub><mi>C</mi><mn>1</mn></msub>")
    end
    print(io, "")
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
            push!(c_strings, "<msubsup><mi>C</mi><mn>$key</mn><mn>$value</mn></msubsup>")
        else
            push!(c_strings, "<msub><mi>C</mi><mn>$key</mn></msub>")
        end
    end
    print(io, join(c_strings, "<mo>&times;</mo>"))
    print(io, " ")
end



preamble = """<!DOCTYPE html>
<html lang="en">
<head>
    <title>SIC data tables</title>
    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1">
    <meta name="keywords" content="SIC, SICPOVM, algebraic number theory, Zauner's conjecture, equiangular lines, Stark conjectures">
</head>
<style>
    td {
        padding-right: 15px;
        padding-top: 3px;
        padding-bottom: 3px;
    }
    th {
        padding: 15px;
        text-align: center;
    }
    thead {
        border-bottom: thin solid #000000;
        border-top: medium solid #000000;
    }
    tbody {
        border-bottom: medium solid #000000;
    }
    tr, td .third {
        text-align: right;
    }
    tr:hover {background-color: #f5f5f5}
    p {text-align: justify;}
    table {
        border-collapse: collapse;
    }
</style>

<body style="margin:25px;width:800px;margin-bottom:100px; ">

<h1>SIC data tables</h1>

<p>This table contains algebraic data for each inequivalent admissible tuple in dimensions <math><mn>4</mn><mo>–</mo><mn>100</mn></math>, comprising <math><mn>614</mn></math> total tuples.
This list is conjecturally complete for all Weyl-Heisenberg covariant SICs in these dimensions.
Each admissible tuple is specified by a dimension <math><mi>d</mi></math> and an integer binary quadratic form <math><mi>Q</mi></math> as follows.
First factorize <math>
  <mrow>
  <mo stretchy="false">(</mo>
  <mi>d</mi>
  <mo>+</mo>
  <mn>1</mn>
  <mo stretchy="false">)</mo>
  <mo stretchy="false">(</mo>
  <mi>d</mi>
  <mo>-</mo>
  <mn>3</mn>
  <mo stretchy="false">)</mo>
  <mo> =</mo>
  <msup>
    <mi>s</mi>
    <mn>2</mn>
  </msup>
  <msub>
    <mi>&Delta;</mi>
    <mn>0</mn>
  </msub>
  </mrow>
</math> where <math>
  <msub>
    <mi>&Delta;</mi>
    <mn>0</mn>
  </msub>
</math> is a fundamental discriminant.
Then <math>
  <mrow>
  <mo stretchy="false">(</mo>
  <mi>d</mi>
  <mo>,</mo>
  <mi>Q</mi>
  <mo stretchy="false">)</mo>
  </mrow>
</math> gives an admissible tuple if <math>
<mrow>
  <mi mathvariant="normal">disc</mi>
  <mo stretchy="false">(</mo>
  <mi>Q</mi>
  <mo stretchy="false">)</mo>
  <mo> =</mo>
  <msup>
    <mi>f</mi>
    <mn>2</mn>
  </msup>
  <msub>
    <mi>&Delta;</mi>
    <mn>0</mn>
  </msub>
</mrow>
</math> where <math><mi>f</mi></math> divides <math><mi>s</mi></math>.
The other columns can be computed from these data, but they may be difficult to compute,
for example requiring integer factoring or finding a fundamental unit.
The column <math><msub><mi>&Delta;</mi><mn>0</mn></msub></math> contains the fundamental discriminant of <math><mi>Q</mi></math> and <math><mi>h</mi></math> is the order of the class group
<math>
  <mrow>
  <mi mathvariant="normal">Cl</mi>
  <mo stretchy="false">(</mo>
  <msub><mi mathvariant="script">O</mi><mi>f</mi></msub>
  <mo stretchy="false">)</mo>
  </mrow>
</math>, given in the next two columns respectively.
The Galois group
<math>
  <mrow>
  <mi mathvariant="normal">Gal</mi>
  <mo stretchy="false">(</mo>
  <mi>E</mi>
  <mo>/</mo>
  <mi>K</mi>
  <mo stretchy="false">)</mo>
  </mrow>
</math>
of the field containing the overlaps over
<math>
  <mrow>
  <mi>K</mi>
  <mo> =</mo>
  <mi mathvariant="double-struck">Q</mi>
  <mo stretchy="false">(</mo>
  <msqrt>
    <msub>
      <mi>&Delta;</mi>
      <mn>0</mn>
    </msub>
  </msqrt>
  <mo stretchy="false">)</mo>
  </mrow>
</math> is given in the next column.
As both the class group and the Galois group are finite and abelian, we give the canonical decomposition into cyclic groups <math><msub><mi>C</mi><mi>k</mi></msub></math> of order <math><mi>k</mi></math>.
For the special case that the tuple has so-called <math><msub><mi>F</mi><mi>a</mi></msub></math> symmetry, we have not yet worked out the Galois groups, so we mark these entries as tbd.
The column <math><msup><mi>L</mi><mi>n</mi></msup></math> contains a generator <math><mi>L</mi></math> of the stability group of <math><mi>Q</mi></math> in
<math>
<mrow>
  <msub>
    <mi mathvariant="normal">GL</mi>
    <mi>2</mi>
  </msub>
  <mo stretchy="false">(</mo>
  <mi mathvariant="double-struck">Z</mi>
  <mo stretchy="false">)</mo>
</mrow>
</math>
and its order <math><mi>n</mi></math> in
<math>
<msub>
  <mi mathvariant="normal">GL</mi>
  <mi>2</mi>
</msub>
<mo stretchy="false">(</mo>
<mi mathvariant="double-struck">Z</mi>
<mo>/</mo>
<mi>d</mi>
<mo>'</mo>
<mi mathvariant="double-struck">Z</mi>
<mo stretchy="false">)</mo>
</math>;
that is, treating <math><mi>Q</mi></math> as a symmetric matrix we have
<math>
<mrow>
  <msup>
    <mi>L</mi>
    <mi>T</mi>
  </msup>
  <mi>Q</mi>
  <mi>L</mi>
  <mo> =</mo>
  <mi>det</mi>
  <mo stretchy="false">(</mo>
  <mi>L</mi>
  <mo stretchy="false">)</mo>
  <mi>Q</mi>
</mrow>
</math>
and
<math>
<mrow>
  <msup>
    <mi>L</mi>
    <mi>n</mi>
  </msup>
  <mo> =</mo>
  <mn>1</mn>
  <mo>&VeryThinSpace;</mo>
  <mo stretchy="false">(</mo>
  <mi>mod</mi>
  <mo>&VeryThinSpace;</mo>
  <mi>d</mi>
  <mo>'</mo>
  <mo stretchy="false">)</mo>
</mrow>
</math>.
If the tuple has antiunitary symmetry, we denote this with a Y in the a.u. column.
Finally, <math><mi>&ell;</mi></math> is the length of the word expansion of <math><msup><mi>L</mi><mi>n</mi></msup></math> using the Hirzebruch-Jung (negative regular) reduction into the standard
(<math><mi>S</mi></math> and <math><mi>T</mi></math>) generators of
<math>
<mrow>
  <msub>
    <mi mathvariant="normal">SL</mi>
    <mn>2</mn>
  </msub>
  <mo stretchy="false">(</mo>
  <mi mathvariant="double-struck">Z</mi>
  <mo stretchy="false">)</mo>
</mrow>
</math>.
This is one measure of the complexity of constructing the actual fiducial vector for that input.
The <math><mi>Q</mi></math> in this list were selected among class representatives to minimize this complexity, although this choice is not generally unique.
The data here are sufficient to compute a ghost fiducial in each class,
but to fully specify a SIC,
one must additionally choose a sign-switching Galois automorphism
<math>
<mrow>
  <msqrt>
    <msub>
      <mi>&Delta;</mi>
      <mn>0</mn>
    </msub>
  </msqrt>
  <mo>&#x2192;</mo>
  <mo>-</mo>
  <msqrt>
    <msub>
      <mi>&Delta;</mi>
      <mn>0</mn>
    </msub>
  </msqrt>
</mrow>

</math>
over an appropriate field extension of <math><mi>K</mi></math>.
</p>

<div id="SICs">
  <table id="table" style="margin-left:100px; text-align:left">
    <thead>
      <tr>
        <th><math><mi>d</mi></math></th>
        <th><math><msub><mi>&Delta;</mi><mn>0</mn></msub></math></th>
        <th><math><mi>f</mi></math></th>
        <th><math><mi>h</mi></math></th>
        <th><math><mrow><mi mathvariant="normal">Cl</mi><mo  stretchy="false">(</mo><msub><mi mathvariant="script">O</mi><mi>f</mi></msub><mo stretchy="false">)</mo></mrow></math></th>
        <th><math><mrow><mi mathvariant="normal">Gal</mi><mo  stretchy="false">(</mo><mi>E</mi><mo>/</mo><mi>K</mi><mo stretchy="false">)</mo></mrow></math></th>
        <th><math><mi>Q</mi></math></th>
        <th><math><mi>L</mi></math></th>
        <th><math><mi>n</mi></math></th>
        <th><math><mi>a.u.</mi></math></th>
        <th><math><mi>&ell;</mi></math></th>
      </tr>
    </thead>
      <tbody>
"""

open(file_name, "w") do io
    print(io, preamble)
    d = 1
    q = 0
    # for j=1:3292 # 4 ≤ d ≤ 256
    # for j=1:282 # 4 ≤ d ≤ 64
    # for j = 1:100 # 4 ≤ d ≤ 35
    # for j=1:35 # testing
    for j=1:614 # 4 ≤ d ≤ 100
        F = AdmissibleTuple(dq(j))
        au, L = is_antiunitary_with_generator(F)
        if !(au)
            L = F.L
        end
        print(io,"""        <tr>\n""")
        print(io,"""          <td data-title="d"><math><mn>$( F.d==d ? "" : Int(F.d))</mn></math></td>\n""")
        print(io,"""          <td data-title="Delta"><math><mn>$( F.d==d ? "" : Int(F.D))</mn></math></td>\n""")
        print(io,"""          <td data-title="f"><math><mn>$( F.d==d && F.q == q ? "" : Int(F.q))</mn></math></td>\n""")
        if !(F.d == d && F.q == q)
            cgp = class_group_structure(F.q^2 * F.D)
            print(io,"""          <td data-title="h"><math><mn>$(prod(cgp))</mn></math></td>\n""") # class number
            print(io,"""          <td data-title="Cl(O_f)"><math>""")
            _print_cycle_decomp(io, cgp)
            print(io,"""</math></td>\n""")
        else
            print(io, """          <td data-title="h"></td>\n""")
            print(io, """          <td data-title="Cl(O_f)"></td>\n""")
        end
        if !(F.d == d && F.q == q)
            print(io,"""          <td data-title="Gal(E/K)">""")
            if has_fa_symmetry(F)
                print(io, "<math><mn> tbd </math></mn>") # galois group
            else
                ords = galois_group_structure(F)
                print(io,"""<math>""")
                _print_cycle_decomp(io, ords) # galois group
                print(io,"""</math>""")
            end
            print(io, """</td>\n""")
        else
            print(io,"""          <td data-title="Gal(E/K)"></td>""")
        end
        print(io,"""          <td data-title="Q">&langle;<math><mn>$(Int(F.Q.a))</mn><mi>,</mi><mn>$(Int(F.Q.b))</mn><mi>,</mi><mn>$(Int(F.Q.c))</mn></math>&rangle;</td>\n""")
        print(io,"""          <td data-title="L"><math><mrow><mo>(</mo>\n""")
        print(io,"""            <mtable rowspacing="4px" columnspacing="6px" columnalign="center">\n""")
        print(io,"""              <mtr><mtd><mn>$(Int(L[1,1]))</mn></mtd><mtd><mn>$(Int(L[1,2]))</mn></mtd></mtr>\n""")
        print(io,"""              <mtr><mtd><mn>$(Int(L[2,1]))</mn></mtd><mtd><mn>$(Int(L[2,2]))</mn></mtd></mtr>\n""")
        print(io,"""            </mtable><mo>)</mo></mrow></math>\n""")
        print(io,"""          </td>\n""")
        print(io,"""          <td data-title="n"><math><mn>$(Int(2^au*F.k))</mn></math></td>\n""")
        print(io,"""          <td data-title="antiunitary"><math><mn>$(au ? "Y" : "")</mn></math></td>\n""")
        print(io,"""          <td data-title="ell"><math><mn>$(Int(length(psl2word(F.A))))</mn></math></td>\n""")
        print(io,"""        </tr>\n""")
        d = F.d
        q = F.q
    end
    epilogue = """
      </tbody>
    </table>
    </div>

    </body>
</html>
"""
    print(io,epilogue)
end
