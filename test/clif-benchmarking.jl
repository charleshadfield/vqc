using VirtualQuantumComputer, Test

n = 10
q1 = rand(1:n)
q2 = rand(1:n)
while q1 == q2
    q2 = rand(1:n)
end

vqc = VQC(n)

symp = buildsymplecticrepn(vqc)

@test symp[q1, q1] == 1
@test symp[q2, q2] == 1


oneQclifford!(symp, :H, q1)
twoQclifford!(symp, :CNOT, q1, q2)

columnindices = [q1, q2, n+q1, n+q2, 2n+1]
reducedsymp = symp[[q1,q2], columnindices]

bell = buildsymplecticrepn(VQC(2))
oneQclifford!(bell, :H, 1)
twoQclifford!(bell, :CNOT, 1, 2)


@test reducedsymp == bell
