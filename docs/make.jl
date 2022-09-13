push!(LOAD_PATH,"../src/maxwell","../examples/vanderPol","../src/wave1D","../src/nbody","../src/TOV")

using Documenter, Maxwell, Wave1D, NBody

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
                    "N-Body Gravity" => "nbody.md",
                    "Wave Eq. 1D" => "wave1D.md",
                    "Maxwell" => "maxwell.md"
                  ]
               ]
)

deploydocs(
    repo = "github.com/davidneilsen/numex.git",
)
