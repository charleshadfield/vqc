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
    @test rref!(sym) == symreduced

    sym = buildsymplecticrepn(VQC(3))
    sym[1, 2] = 1
    sym[2, 7] = 1
    symreduced = buildsymplecticrepn(VQC(3))
    symreduced[1, 7] = 1
    symreduced[2, 7] = 1
    @test rref!(sym) == symreduced
end
