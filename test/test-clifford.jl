
@testset "symplectic construction" begin
    vqc1 = VQC(1)
    sym1 = buildsymplecticrepn(vqc1)
    @test sym1 == [1 0 0]

    vqc3 = VQC(3)
    sym3 = buildsymplecticrepn(vqc3)
    @test sum(sym3) == 3
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

@testset "reduced twoQclifford!" begin
    sym = [0; 1; 0; 0] # IZ
    ZZ = [1; 1; 0; 0]
    twoQclifford!(sym, :CNOT)
    @test sym == ZZ

    #sym = [1; 0; 0; 1] #ZX
    #YY = [1; 1; 1; 1]
    #@test twoQclifford!(sym, :CNOTreversed) == YY

    sym = [0; 0; 1; 1]
    II = [0; 0; 0; 0]
    @test twoQclifford!(sym, :CZ) == II
end

@testset "twoQclifford!" begin
    sym1 = buildsymplecticrepn(VQC(2))
    twoQclifford!(sym1, :CNOT, 1, 2)
    sym2 = buildsymplecticrepn(VQC(2))
    sym2[2, 1] = 1
    @test sym1 == sym2

    sym1 = buildsymplecticrepn(VQC(2))
    twoQclifford!(sym1, :CNOT, 2, 1)
    sym2 = buildsymplecticrepn(VQC(2))
    sym2[1, 2] = 1
    @test sym1 == sym2

    bell = [1 1 0 0 0; 0 0 1 1 0]
    sym = buildsymplecticrepn(VQC(2))
    oneQclifford!(sym, :H, 1)
    twoQclifford!(sym, :CNOT, 1, 2)
    VirtualQuantumComputer.rref!(sym)
    @test sym == bell
end
