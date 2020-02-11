push!(LOAD_PATH,"../src/")
using Documenter, MolHandler
makedocs(
    sitename="MolHandler.jl",
	  pages = [
        "Top" => "index.md",
        "Functions" => "functions.md",
        "Structs" => "structs.md"
    ]
)
