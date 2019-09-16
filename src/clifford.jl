export buildsymplecticrepn, oneQclifford!, twoQclifford!

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

"""
    twoQgate!(reducedsymp, gate::Symbol)

Apply `CNOT` or `CZ` to `reducedsymp` where
`reducedsymp = [z1; z2; x1; x2;]`. Here `zi`, `xi` are binary, and `reducedsymp` should
be considered a **column** vector. So multiplication is from the right.

Note that both `CZ` and `CNOT` do not required knowledge of signed bit.
Moreover, they are linear in the sense that we only need to define their
actions on the `z`-part and `x`-part independently.
"""
function twoQclifford!(reducedsymp, gate::Symbol)
    @assert length(reducedsymp) == 4
    @assert gate in [:CNOT, :CZ]
    if gate == :CNOT
        reducedsymp[1:2] = clifCNOTpartZ * reducedsymp[1:2]
        reducedsymp[3:4] = clifCNOTpartX * reducedsymp[3:4]
    else
        reducedsymp[3:4] = clifCZpartX * reducedsymp[3:4]
    end
    # reduce modulo 2
    reducedsymp[:] = reducedsymp .% 2
    reducedsymp
end

"""
    twoQclifford!(symp, gate::Symbol, q1::Int, q2::Int)

Apply clifford `gate` to `symp` by repeatedly calling
`twoQclifford!(reducedsymp, gate)`
`reducedsymp = [z1; z2; x1; x2]` for binary values `zi`, `xi`.
"""
function twoQclifford!(symp, gate::Symbol, q1::Int, q2::Int)
    n = size(symp)[1]
    @assert q1 ≤ n
    @assert q2 ≤ n
    @assert gate in [:CNOT, :CZ]

    # the following code works for CNOT even if q1 > q2
    for row = 1:n
        reducedsymp = symp[row, [q1, q2, n+q1, n+q2]]
        symp[row, [q1, q2, n+q1, n+q2]] = twoQclifford!(reducedsymp, gate)
    end
    symp
end
