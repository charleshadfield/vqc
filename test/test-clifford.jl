using VirtualQuantumComputer: rref!

@testset "symplectic construction" begin
    vqc1 = VQC(1)
    sym1 = buildsymplecticrepn(vqc1)
    @test sym1 == [1 0 0]

    vqc3 = VQC(3)
    sym3 = buildsymplecticrepn(vqc3)
    @test sum(sym3) == 3
end

@testset "rref!" begin
    sym = [1 1 0; 0 1 0]
    symreduced = [1 0 0; 0 1 0]
    rref!(sym)
    @test sym == symreduced

    sym = buildsymplecticrepn(VQC(3))
    sym[1, 2] = 1
    sym[2, 7] = 1
    symreduced = buildsymplecticrepn(VQC(3))
    symreduced[1, 7] = 1
    symreduced[2, 7] = 1
    @test rref!(sym) == symreduced
end

@testset "reduced oneQclifford!" begin
    sym = buildsymplecticrepn(VQC(1))
    oneQclifford!(sym, :X)
    @test sym == [1 0 1]

    @test oneQclifford!(sym, :H) == [0 1 1]
    @test oneQclifford!(sym, :Z) == [0 1 0]
    @test oneQclifford!(sym, :S) == [1 1 0]
    @test oneQclifford!(sym, :Y) == [1 1 0]
    @test oneQclifford!(sym, :Z) == [1 1 1]
    oneQclifford!(sym, :S)
    oneQclifford!(sym, :H)
    @test sym == [1 0 0]
end

@testset "oneQclifford!" begin
    sym = buildsymplecticrepn(VQC(1))
    oneQclifford!(sym, :X, 1)
    @test sym == [1 0 1]

    sym = buildsymplecticrepn(VQC(2))
    oneQclifford!(sym, :X, 1)
    oneQclifford!(sym, :H, 2)
    outcome = [[1 0 0 0 1]; [0 0 0 1 0]]
    @test sym == outcome
end


#bell = [[1 1 0 0 0]; [0 0 1 1 0]]
