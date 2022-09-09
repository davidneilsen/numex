push!(LOAD_PATH,"../src/maxwell","../examples/vanderPol")

using Documenter, Maxwell

makedocs(
         sitename="numex",
         format = Documenter.HTML(
             # Use clean URLs, unless built as a "local" build
             # prettyurls = !("local" in ARGS)
             prettyurls = false
         ),
         pages=[
                 "Home" => "index.md",
                 "Examples" => Any[
                    "Nonlinear Oscillator" => "vanderpol.md",
                    "Maxwell" => "maxwell.md"
                  ]
               ]
)

deploydocs(
    repo = "github.com/davidneilsen/numex.git",
)
