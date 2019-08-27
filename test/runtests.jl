using VirtualQuantumComputer, Test

@testset "All tests" begin

    @testset "trivial" begin
        @test true
    end

    include("test-gates.jl")

    include("test-measurement.jl")

    include("test-clifford.jl")

    include("test-clifford-measurement.jl")

end # All tests.
