using VirtualQuantumComputer, Test

@testset "trivial" begin
    @test true
end

include("test-gates.jl")

include("test-measurement.jl")

include("test-clifford.jl")
