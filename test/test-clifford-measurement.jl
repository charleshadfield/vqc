using VirtualQuantumComputer: rref!

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

    #bell
    sym         = [1 1 1 1 0; 0 0 1 1 0]
    symreduced  = [1 1 0 0 1; 0 0 1 1 0]
    @test rref!(sym) == symreduced
end

@testset "measureclifford! simple" begin
    vqc = VQC(1)
    symp = buildsymplecticrepn(vqc)

    measurement = measureclifford!(symp, "X", [1])
    symp
    @test symp[1:2] == [0, 1] && symp[3] == measurement

    vqc = VQC(1)
    symp = buildsymplecticrepn(vqc)
    oneQclifford!(symp, :X, 1)
    measurement = measureclifford!(symp, "Z", [1])
    @test symp == [1 0 1] && measurement == 1

    vqc = VQC(2)
    symp = buildsymplecticrepn(vqc)
    oneQclifford!(symp, :X, 1)
    measurement = measureclifford!(symp, "ZZ", [1, 2])
    @test measurement == 1 && symp == [1 0 0 0 1; 0 1 0 0 0]

    vqc = VQC(2)
    symp = buildsymplecticrepn(vqc)
    symp[1,2] = 1
    measurement = measureclifford!(symp, "ZZ", [1, 2])
    @test measurement == 0 && symp == [1 0 0 0 0; 0 1 0 0 0]
end

@testset "measureclifford! medium" begin

    # bell with Y measurement
    sym = buildsymplecticrepn(VQC(2))
    oneQclifford!(sym, :H, 1)
    twoQclifford!(sym, :CNOT, 1, 2)
    @test sym == [0 0 1 1 0; 1 1 0 0 0]
    m1 = measureclifford!(sym, "Y", [1])
    sym
    m2 = measureclifford!(sym, "Y", [2])
    @test sym == [1 0 1 0 m1; 0 1 0 1 m2]

    # alternate bell with YY measurement
    sym = buildsymplecticrepn(VQC(2))
    oneQclifford!(sym, :H, 1)
    twoQclifford!(sym, :CNOT, 1, 2)
    @test sym == [0 0 1 1 0; 1 1 0 0 0]

    measurement = measureclifford!(sym, "YY", [1, 2])
    @test measurement == 1
end
