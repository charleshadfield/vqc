export VQC, buildwf, oneQgate!, twoQgate!, measure!

struct VQC
    n::Int
end

function buildwf(vqc::VQC)
    # 2*2*2*...*2 Array
    wf = zeros(ComplexF64, (2 for _ in 1:vqc.n)...)
    # initialize to |0, 0, ..., 0> state
    wf[ones(Int, vqc.n)...] = 1
    wf
end

function oneQgate(reducedwf, A)
    # subroutine of A acting on single-qubit-wave-function
    A * reducedwf
end

function oneQgate!(wf, A, q::Int)
    n = ndims(wf)
    @assert q ≤ n
    for dummyindex in Iterators.product((1:2 for _ in 1:n-1)...)
        a = dummyindex[1:q-1]
        b = dummyindex[q:n-1]
        wf[a..., :, b...] = oneQgate(wf[a..., :, b...], A)
    end
    wf
end

function twoQgate(reducedwf, A, bool::Bool)
    # subroutine of A acting on two-qubit-wave-function
    # :bool: T if in order, F if out of order
    newreducedwf = zeros(ComplexF32, 2, 2)
    if bool
        for (i, j) in Iterators.product((1:2), (1:2))
            newreducedwf[i, j] = sum(A[i,j,x,y]*reducedwf[x,y] for (x,y) in Iterators.product((1:2), (1:2)))
        end
    else
        for (i, j) in Iterators.product((1:2), (1:2))
            newreducedwf[i, j] = sum(A[j,i,x,y]*reducedwf[y,x] for (x,y) in Iterators.product((1:2), (1:2)))
        end
    end
    newreducedwf
end

function twoQgate!(wf, A, q1::Int, q2::Int)#::Array{Complex{Float}, 2}, wf::Array{Complex{Float}})
    n = ndims(wf)
    @assert n ≥ 2
    @assert q1 != q2
    @assert q1 ≤ n
    @assert q2 ≤ n
    qmin = min(q1, q2)
    qmax = max(q1, q2)
    q1 < q2 ? bool = true : bool = false
    for dummyindex in Iterators.product((1:2 for _ in 1:n-2)...)
        a = dummyindex[1:qmin-1]
        b = dummyindex[qmin:qmax-2]
        c = dummyindex[qmax-1: n-2]
        wf[a..., :, b..., :, c...] = twoQgate(wf[a..., :, b..., :, c...], A, bool)
    end
    wf
end

function probabilities(wf, q::Int)
    # obtain two probabilities associated with measuring 0 or 1 in Z basis
    n = ndims(wf)
    @assert q ≤ n
    a = repeat([:], q-1)
    b = repeat([:], n-q)
    wf[a..., 1, b...]
    prob0 = sum(abs2(x) for x in wf[a..., 1, b...])
    prob1 = sum(abs2(x) for x in wf[a..., 2, b...])
    #@assert prob0 + prob1 ≈ 1
    #understand complex/real trade off when asserting approx equality
    [prob0, prob1]
end

function collapse!(wf, q::Int, measurement::Int, prob)
    # rescale wf by sqrt prob in appropriate index
    # set to zero entries in otherindex
    # measurement is either 0 or 1
    # return modified wf
    index = measurement + 1
    otherindex = (measurement + 1) % 2 + 1
    @assert index + otherindex == 3
    n = ndims(wf)
    a = repeat([:], q-1)
    b = repeat([:], n-q)
    wf[a..., index, b...] = 1/√prob * wf[a..., index, b...]
    wf[a..., otherindex, b...] = 0 * wf[a..., otherindex, b...]
    wf
end

function measure!(wf, q::Int)
    prob0, prob1 = probabilities(wf, q)
    random = rand()
    if random < prob0
        collapse!(wf, q, 0, prob0)
        0
    else
        register = 1
        collapse!(wf, q, 1, prob1)
        1
    end
end
