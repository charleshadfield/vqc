println("Hello")
cd("/Users/charles/Desktop/julia/vqc")
println(pwd())
using Pkg
Pkg.activate(".")

include("base.jl")

using Test

include("test-gates.jl")

include("test-measurement.jl")
