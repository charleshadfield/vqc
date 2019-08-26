using VirtualQuantumComputer, Test

function measureclifford!(symp, pauli::String, q::Int)
end

function measureclifford!(symp, pauli::String, qubits::Vector{Int})
    n = size(symp)[1]
    @assert length(string) == length(qubits)
    weight = length(string)
    for k = 1:weight
        @assert qubits[k] ≤ n
        @assert string[k] in ["X", "Y", "Z"]
    end

    # if qubits are out of order, then reorder
    if sort(qubits) != qubits
        permutation = sortperm(qubits)
        measureclifford!(symp, string[permutation], qubits[permutation])
    end

    # build 2*weight length vector representing measurement in symplectic form
    measurementrepn = zeros(Int, 2*weight)
    for k = 1:weight
        pauli[k] == "X" && measurementrepn[n+k] = 1
        pauli[k] == "Y" && measurementrepn[k] = measurementrepn[n+k]= 1
        pauli[k] == "Z" && measurementrepn[k] = 1
    end

    # build 2*weight length vector to see which stabilizers anticommute
    # with pauli measurement
    # this is the measurementrepn upon application of symplectic form
    dualmeasurementrepn = zeros(Int, 2*weight)
    for k = 1:weight
        pauli[k] == "X" && measurementrepn[k] = 1
        pauli[k] == "Y" && measurementrepn[k] = measurementrepn[n+k]= 1
        pauli[k] == "Z" && measurementrepn[n+k] = 1
    end

    # find noncommuting stabilizers
    noncommutingrows::Vector{Int} = []
    paulipositions = vcat(qubits, n .+ qubits)
    for row = 1:n
        parity = (symp[row, paulipositions] .* dualmeasurementrepn) .% 2
        parity == 1 && append!(noncommutingrows, row)
    end

    # if non-commuting rows exist then measurement is random
    # set random variable measurementoutcome

    # must update non-commuting rows by multiplying pairs of rows together.
    # e.g. if rows 1,2,3 have stabilizers S1,S2,S3 anti-commuting with
    #      pauli measurement, then we perform the following
    #      rows 2,3 are updated to be S1*S2, S1*S3
    #      row 1 is replaced by paulimeasurement and
    #      the parity of row 1 is measurmentoutcome
    if length(noncommutingrows) != 0
        measurementoutcome = rand(0:1)
        # update noncommuting rows S2, ...
        for row in noncommutingrows[2:]
            symp[row, :] = (symp[row, :] .* symp[1, :]) .% 2
        end
        # update noncommuting row S1 with appropriately chosen measurement
        newstabilizer = zeros(Int, 2n+1)
        newstabilizer[paulipositions] = measurementrepn
        newstabilizer[end] = measurementoutcome
        symp[noncommutingrows[1]] = newstab

    # if no non-commuting rows exist then measurement is determinisitc.
    #
    else
        println("measurement is deterministic")
    end
    println("not yet finished")
    measurementoutcome
end
