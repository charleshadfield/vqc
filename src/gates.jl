export X, Y, Z, H, S, CNOT, CZ

X = [0 1; 1 0]
Y = [0 -im; im 0]
Z = [1 0; 0 -1]

H = 1/√2 * (X + Z)
S = [1 0; 0 im]


CNOT = zeros(Int, 2, 2, 2, 2)
CNOT[1,1,1,1] = 1
CNOT[1,2,1,2] = 1
CNOT[2,2,2,1] = 1
CNOT[2,1,2,2] = 1
CNOT

CZ = zeros(Int, 2, 2, 2, 2)
CZ[1,1,1,1] = 1
CZ[1,2,1,2] = 1
CZ[2,1,2,1] = 1
CZ[2,2,2,2] = -1
CZ
