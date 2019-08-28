export measureclifford!


"""
    paulimultiplication!(paulis::Matrix{Int})

`paulis = [z1 x1; z2 x2]` represents two paulis ``P_1, P_2``.
Function modifies paulisquare to be
`[z3 x3; z2 x2]`
where ``P_3 = P_1 \\circ P_2`` and returns `phase` ``\\in \\{-1,0,1\\}``.

e.g.
- ``P_1 = X, P_2 = X``
    - original `paulis = [0 1; 0 1]`
    - modified `paulis = [0 1; 0 1]`
    - return `phase = 0`
- ``P_1 = X, P_2 = Y``
    - original `paulis = [0 1; 1 1]`
    - modified `paulis = [1 0; 1 1]`
    - return `phase = 1` since ``X\\circ Y = i Z = (i)^1 Z``
"""
function paulimultiplication!(paulis::Matrix{Int})
    phase = 0
    det = (paulis[1,1]*paulis[2,2] + paulis[1,2]*paulis[2,1]) % 2
    if det != 0
        paulis == [0 1; 1 1] && (phase = 1)
        paulis == [0 1; 1 0] && (phase = -1)
        paulis == [1 1; 0 1] && (phase = -1)
        paulis == [1 1; 1 0] && (phase = 1)
        paulis == [1 0; 0 1] && (phase = 1)
        paulis == [1 0; 1 1] && (phase = -1)
    end
    paulis[1,:] = (paulis[1,:] .+ paulis[2,:]) .% 2
    phase
end



"""
    phasedstabilizermultiplication!(stabs, n::Int)

- `stabs = [stab1; stab2]` represents two stabilizers.
- `n` is number of qubits

`stab1` updated to be `stab1 * stab2` with appropriate phase.
Function repeatedly calls function `paulimultiplication!`
"""
function phasedstabilizermultiplication!(stabs::Matrix{Int}, n::Int)
    @assert size(stabs) == (2, 2n+1)
    phasemod4 = 0
    for k = 1:n
        paulis = stabs[:,[k, n+k]]
        phasemod4 += paulimultiplication!(paulis)
        stabs[:,[k, n+k]] = paulis
    end
    @assert phasemod4 % 4 in [-2, 0, 2]
    sign = Int(abs(phasemod4 % 4) // 2)
    stabs[1,2n+1] = (stabs[1,2n+1] + stabs[2,2n+1] + sign) % 2
    stabs
end

"""
    rref!(symp)

Row reduce symplectic matrix.

Algorithm assumes entries are binary valued.
"""
function rref!(symp)
    # get dimensions of matrix, final column is pm1 value
    nr, nc = size(symp) .+ (0, -1)
    @assert 2nr == nc || 2(nr-1) == nc
    # If 2nr = nc then rref! is being used on stabilizer generators as usual.
    # Otherwise an extra row has been added in order to calculate measurement.
    2nr == nc ? (nqubits = nr) : (nqubits = nr - 1)
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
                stabs = symp[[iprime, i], :]
                phasedstabilizermultiplication!(stabs, nqubits)
                symp[[iprime, i], :] = stabs
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
    @assert length(pauli) == length(qubits)
    weight = length(pauli)
    for k = 1:weight
        @assert qubits[k] ≤ n
        @assert pauli[k] in ['X', 'Y', 'Z']
    end

    # if qubits are out of order, then reorder
    if sort(qubits) != qubits
        permutation = sortperm(qubits)
        measureclifford!(symp, string[permutation], qubits[permutation])
    end

    # build 2*weight length vector representing measurement in symplectic form
    measurementrepn = zeros(Int, 2*weight)
    for k = 1:weight
        pauli[k] == 'X' && (measurementrepn[weight+k] = 1)
        pauli[k] == 'Y' && (measurementrepn[k] = measurementrepn[weight+k]= 1)
        pauli[k] == 'Z' && (measurementrepn[k] = 1)
    end

    # build 2*weight length vector to see which stabilizers anticommute
    # with pauli measurement
    # this is the measurementrepn upon application of symplectic form
    dualmeasurementrepn = zeros(Int, 2*weight)
    for k = 1:weight
        pauli[k] == 'X' && (dualmeasurementrepn[k] = 1)
        pauli[k] == 'Y' && (dualmeasurementrepn[k] = dualmeasurementrepn[weight+k]= 1)
        pauli[k] == 'Z' && (dualmeasurementrepn[weight+k] = 1)
    end

    # find noncommuting stabilizers
    noncommutingrows::Vector{Int} = []
    paulipositions = vcat(qubits, n .+ qubits)
    for row = 1:n
        parity = sum(symp[row, paulipositions] .* dualmeasurementrepn) % 2
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
            symp[row, :] = (symp[row, :] .+ symp[1, :]) .% 2
        end
        # update noncommuting row S1 with appropriately chosen measurement
        newstabilizer = zeros(Int, 2n+1)
        newstabilizer[paulipositions] = measurementrepn
        newstabilizer[end] = measurementoutcome
        symp[noncommutingrows[1],:] = newstabilizer

    # if no non-commuting rows exist then measurement is determinisitc.
    # can find the measurementoutcome by the following trick:
    #   - augment the symplectic representation with a final line which is the
    #     pauli measurement
    #   - initialize the parity of the final stabilizer to "0"
    #   - perform rref! on the augmented symplectic representation
    else
        augmentedsymp = zeros(Int64, n+1, 2n+1)
        augmentedsymp[1:end-1, :] = symp
        augmentedsymp[end, paulipositions] = measurementrepn
        rref!(augmentedsymp)
        measurementoutcome = augmentedsymp[end,end]
        # maybe set symplectic representation to rref! version of original symplectic repn.
        symp[:,:] = augmentedsymp[1:end-1, :]
    end
    measurementoutcome
end
