




"""
    rref!(symp)

Row reduce symplectic matrix.

Algorithm assumes entries are binary valued.
"""
function rref!(symp)
    # get dimensions of matrix, final column is pm1 value
    nr, nc = size(symp) .+ (0, -1)
    i = j = 1
    while i <= nr && j <= nc
        entry, pos = findmax(symp[i:nr, j])
        if entry == 0
            # no leading 1 will be found in column j
            # there may be 1s appearing above row i in column j
            j += 1
            continue
        end
        # a leading 1 has been found at position (i+pos-1, j)
        if pos != 1
            symp[i, :], symp[i+pos-1, :] = symp[i+pos-1, :], symp[i, :]
        end
        # (i, j) is now a leading entry
        # remove all appearances of 1 in column j except in (i, j)
        # do this by adding mod 2 row i to row iprime
        for iprime in 1:nr
            if iprime != i && symp[iprime, j] == 1
                for jprime in j:nc+1
                    symp[iprime, jprime] =  (symp[iprime, jprime] + symp[i, jprime]) % 2
                end
            end
        end
        i += 1
        j += 1
    end
    symp
end




#function measureclifford!(symp, pauli::String, q::Int)
#end

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
        pauli[k] == "X" && (measurementrepn[n+k] = 1)
        pauli[k] == "Y" && (measurementrepn[k] = measurementrepn[n+k]= 1)
        pauli[k] == "Z" && (measurementrepn[k] = 1)
    end

    # build 2*weight length vector to see which stabilizers anticommute
    # with pauli measurement
    # this is the measurementrepn upon application of symplectic form
    dualmeasurementrepn = zeros(Int, 2*weight)
    for k = 1:weight
        pauli[k] == "X" && (measurementrepn[k] = 1)
        pauli[k] == "Y" && (measurementrepn[k] = measurementrepn[n+k]= 1)
        pauli[k] == "Z" && (measurementrepn[n+k] = 1)
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
        for row in noncommutingrows[2:end]
            symp[row, :] = (symp[row, :] .* symp[1, :]) .% 2
        end
        # update noncommuting row S1 with appropriately chosen measurement
        newstabilizer = zeros(Int, 2n+1)
        newstabilizer[paulipositions] = measurementrepn
        newstabilizer[end] = measurementoutcome
        symp[noncommutingrows[1]] = newstab

    # if no non-commuting rows exist then measurement is determinisitc.
    # can find the measurementoutcome by the following trick:
    #   - augment the symplectic representation with a final line which is the
    #     pauli measurement
    #   - initialize the parity of the final stabilizer to "0"
    #   - perform rref! on the augmented symplectic representation
    else
        println("measurement is deterministic")
    end
    println("not yet finished")
    measurementoutcome
end
