export buildsymplecticrepn, oneQclifford!

include("clif-gates.jl")


"""
    buildsymplecticrepn(vqc::VQC)

Construct an matrix with entries `0,1` which is a binary symplectic representation
of the quantum state using the stabilizer formalism.

Let ``n`` be the number of qubits. The matrix is of dimension ``n \\times 2n+1``
and each row is the representation of a stabilizer ``S``.

A stabilizer ``S`` is a Pauli tensor acting on the quantum state ``\\psi`` such
that ``S\\psi = \\psi``. By declaring ``n`` independent stabilizers we have
uniquely defined the quantum state (up to an unimportant ``\\mathbb{C}^*``
projective factor).

The representation of ``S`` as a ``2n+1`` vector `V` is as follows:
- for ``1\\le k \\le n``:
    - presence of `V[k]` indicates the Pauli ``Z_k``;
    - presence of `V[n+k]` indicates the Pauli ``X_k``;
    - presence of both `V[k]` and `V[n+k]` indicates the Pauli ``Y_k``;
- `V[2n+1]` indicates a global ``\\pm 1`` factor where the integer ``s`` provides
    the factor ``(-1)^s``.
"""
function buildsymplecticrepn(vqc::VQC)
    # rectangular array
    symp = zeros(Int64, vqc.n, 2*vqc.n+1)
    # initialize to |0, 0, ..., 0> state
    for q in 1:vqc.n
        symp[q, q] = 1
    end
    symp
end

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


"""
    oneQclifford!(reducedsymp, gate::Symbol)

Apply clifford `gate` to `reducedsymp` where
`reducedsymp = [z, x, s]` for binary values `z, x, s`.

This is a subroutine on a one qubit system.
"""
function oneQclifford!(reducedsymp, gate::Symbol)
    @assert gate in [:X, :Y, :Z, :H, :S]
    # apply matrix with integer outcomes
    if gate == :X
        reducedsymp[:] = reducedsymp * clifX
    elseif gate == :Y
        reducedsymp[:] = reducedsymp * clifY
    elseif gate == :Z
        reducedsymp[:] = reducedsymp * clifZ
    elseif gate == :H
        zx = reducedsymp[1] * reducedsymp[2]
        reducedsymp[:] = reducedsymp * clifHpart
        reducedsymp[3] += zx
    elseif gate == :S
        zx = reducedsymp[1] * reducedsymp[2]
        reducedsymp[:] = reducedsymp * clifSpart
        reducedsymp[3] += zx
    end
    # reduce modulo 2
    reducedsymp[:] = reducedsymp .% 2
    reducedsymp
end

"""
    oneQclifford!(symp, gate::Symbol, q::Int)

Apply clifford `gate` to `symp` by repeatedly calling
`oneQclifford!(reducedsymp, gate)`
`reducedsymp = [z, x, s]` for binary values `z, x, s`.

"""
function oneQclifford!(symp, gate::Symbol, q::Int)
    n = size(symp)[1]
    @assert q ≤ n
    for row = 1:n
        reducedsymp = transpose(symp[row, [q, n+q, 2n+1]])
        symp[row, [q, n+q, 2n+1]] = oneQclifford!(reducedsymp, gate)
    end
    symp
end
