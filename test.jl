print("Hello")
cd("/Users/charles/Desktop/julia/vqc")
using Pkg
Pkg.activate(".")

include("base.jl")

using Test

@testset "oneQbitflip" begin
    vqc = VQC(2)

    wf00 = buildwf(vqc)
    wf10 = buildwf(vqc)
    wf10[1, 1] = 0
    wf10[2, 1] = 1
    # check that ground goes to excited
    @test oneQgate!(wf00, X, 1) == wf10

    # check that excited returns to ground
    @test oneQgate!(wf10, X, 1) == buildwf(vqc)

    # check that wf00 was appropriately modified by oneQgate!
    @test oneQgate!(wf00, X, 1) == buildwf(vqc)
end

@testset "twoQCNOTsubroutine" begin
    vqc = VQC(2)

    wf00 = buildwf(vqc)
    wf10 = oneQgate!(buildwf(vqc), X, 1)
    wf01 = oneQgate!(buildwf(vqc), X, 2)
    wf11 = oneQgate!(oneQgate!(buildwf(vqc), X, 1), X, 2)

    @test twoQgate(wf00, CNOT, true) == buildwf(vqc)
    @test twoQgate(wf10, CNOT, true) == wf11
    @test twoQgate(wf01, CNOT, true) == oneQgate!(buildwf(vqc), X, 2)
    @test twoQgate(wf11, CNOT, true) == oneQgate!(buildwf(vqc), X, 1)

    @test twoQgate(wf00, CNOT, false) == buildwf(vqc)
    @test twoQgate(wf10, CNOT, false) == oneQgate!(buildwf(vqc), X, 1)
    @test twoQgate(wf01, CNOT, false) == oneQgate!(oneQgate!(buildwf(vqc), X, 1), X, 2)
    @test twoQgate(wf11, CNOT, false) == oneQgate!(buildwf(vqc), X, 2)
end

@testset "twoQCNOTroutine" begin
    vqc = VQC(2)

    wf00 = buildwf(vqc)
    wf10 = oneQgate!(buildwf(vqc), X, 1)
    wf01 = oneQgate!(buildwf(vqc), X, 2)
    wf11 = oneQgate!(oneQgate!(buildwf(vqc), X, 1), X, 2)

    @test twoQgate!(wf00, CNOT, 1, 2) == buildwf(vqc)
    @test twoQgate!(wf10, CNOT, 1, 2) == wf11
    @test twoQgate!(wf01, CNOT, 1, 2) == oneQgate!(buildwf(vqc), X, 2)
    @test twoQgate!(wf11, CNOT, 1, 2) == oneQgate!(buildwf(vqc), X, 1)

    wf00 = buildwf(vqc)
    wf10 = oneQgate!(buildwf(vqc), X, 1)
    wf01 = oneQgate!(buildwf(vqc), X, 2)
    wf11 = oneQgate!(oneQgate!(buildwf(vqc), X, 1), X, 2)

    @test twoQgate!(wf00, CNOT, 2, 1) == buildwf(vqc)
    @test twoQgate!(wf10, CNOT, 2, 1) == oneQgate!(buildwf(vqc), X, 1)
    @test twoQgate!(wf01, CNOT, 2, 1) == oneQgate!(oneQgate!(buildwf(vqc), X, 1), X, 2)
    @test twoQgate!(wf11, CNOT, 2, 1) == oneQgate!(buildwf(vqc), X, 2)
end

@testset "threeQCNOTroutine" begin
    vqc = VQC(3)

    wf100 = oneQgate!(buildwf(vqc), X, 1)
    wf101 = oneQgate!(oneQgate!(buildwf(vqc), X, 1), X, 3)
    @test twoQgate!(wf100, CNOT, 1, 3) == wf101

    wf011 = oneQgate!(oneQgate!(buildwf(vqc), X, 2), X, 3)
    wf001 = oneQgate!(buildwf(vqc), X, 3)
    @test twoQgate!(wf011, CNOT, 3, 2) == wf001
end
