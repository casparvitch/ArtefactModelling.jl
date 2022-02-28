using ArtefactModelling
using Documenter

DocMeta.setdocmeta!(
    ArtefactModelling,
    :DocTestSetup,
    :(using ArtefactModelling);
    recursive = true,
)

makedocs(;
    modules = [ArtefactModelling],
    authors = "Sam Scholten <samcaspar@gmail.com> and contributors",
    repo = "https://github.com/casparvitch/ArtefactModelling.jl/blob/{commit}{path}#{line}",
    sitename = "ArtefactModelling.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://casparvitch.github.io/ArtefactModelling.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(;
    repo = "github.com/casparvitch/ArtefactModelling.jl",
    devbranch = "main",
)
