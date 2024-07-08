arrow[v_] := {Arrowheads[Medium], Arrow[Tube[{{0, 0, 0}, v/Sqrt[3]}]]};
v = {{ 1, 1, 1},{ 1,-1,-1},{-1, 1,-1},{-1,-1, 1}}
t = arrow/@v;
s = {Black, Opacity[.3], Specularity[White, 5], Sphere[{0, 0, 0}]};
p = Graphics3D[ Join[t,{s}], Boxed -> False, Lighting -> DirectionalLight[RGBColor[1, 1, 1], {{5, 5, 4}, {5, 5, 0}}], Background -> None, ImageMargins -> 0, ViewPoint -> {Pi, Pi/4, -1/2}];
Export[DirectoryName[$InputFileName] <> "/assets/logo.png", p, "PNG", Background -> None];
