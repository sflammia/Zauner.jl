push!(LOAD_PATH, "../src")

using Documenter, Zauner

DocMeta.setdocmeta!(
    Zauner,
    :DocTestSetup,
    :(using Zauner);
    recursive=true
)

# repo_url = "https://sflammia.github.io/Zauner.jl"
repo_url = "https://github.com/sflammia/Zauner.jl"

makedocs(
    modules=[Zauner],
    doctest=true,
    # warnonly = :doctest,
    authors="Marcus Appleby, Steven T Flammia, Gene S Kopp",
    repo=repo_url * "/{commit}{path}#{line}",
    sitename="Zauner",
    format=Documenter.HTML(;
        canonical=repo_url,
        repolink=repo_url,
        edit_link=repo_url * "{commit}{path}#{line}",
    ),
    checkdocs=:missing,
    pages=[
        "index.md",
        "List of Functions" => "functions.md",
        "Tables" => "tables.md",
    ],
)
