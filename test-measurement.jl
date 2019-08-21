@testset "probabilities" begin

    for n in 1:3
        for q in 1:n
            wf = buildwf(VQC(n))
            @test probabilities(wf, q) == [1, 0]

            oneQgate!(wf, X, q)
            @test probabilities(wf, q) == [0, 1]
        end
    end

    wf = buildwf(VQC(3))
    oneQgate!(wf, H, 2)
    @test probabilities(wf, 1) ≈ [1, 0]
    @test probabilities(wf, 2) ≈ [0.5, 0.5]

end

@testset "1Qcollapse" begin
    wf = buildwf(VQC(1))
    @test collapse!(wf, 1, 0, 1) == buildwf(VQC(1))
    oneQgate!(wf, H, 1)
    @test collapse!(wf, 1, 0, 0.5) ≈ buildwf(VQC(1))

    oneQgate!(wf, H, 1)
    @test collapse!(wf, 1, 1, 0.5) ≈ oneQgate!(buildwf(VQC(1)), X, 1)
end

@testset "2Qcollapse" begin
    wf1 = buildwf(VQC(2))
    wf2 = buildwf(VQC(2))
    oneQgate!(wf1, H, 2)
    oneQgate!(wf2, H, 2)

    @test collapse!(wf1, 1, 0, 1) == wf2

    wf = buildwf(VQC(2))
    oneQgate!(wf, H, 1)
    @test collapse!(wf, 1, 1, 0.5) ≈ oneQgate!(buildwf(VQC(2)), X, 1)
end

@testset "1Qmeasurement" begin
    wf = buildwf(VQC(1))
    @test 0 == measure!(wf, 1)
    oneQgate!(wf, X, 1)
    @test 1 == measure!(wf, 1)

    function measureafterhadamard()
        vqc = VQC(1)
        wf = buildwf(vqc)
        oneQgate!(wf, H, 1)
        measurement = measure!(wf, 1)
        (measurement, wf)
    end

    for _ in 1:10
        measurement, wf = measureafterhadamard()
        if measurement == 0
            @assert wf ≈ buildwf(VQC(1))
        else
            @assert wf ≈ oneQgate!(buildwf(VQC(1)), X, 1)
        end
    end
    @test true

    function average(n)
        sum(measureafterhadamard()[1] for _ in 1:n) / n
    end

    @test 0.4 < average(10000) < 0.6
end

@testset "multipleQmeasurement" begin
    wf = buildwf(VQC(3))
    oneQgate!(wf, H, 2)
    ro1 = measure!(wf, 1)
    ro3 = measure!(wf, 3)

    @test [ro1, ro3] == [0, 0]

    twoQgate!(wf, CNOT, 2, 3)
    ro2 = measure!(wf, 2)
    ro3 = measure!(wf, 3)
    @test ro2 == ro3
end
