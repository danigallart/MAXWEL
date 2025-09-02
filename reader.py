import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import bicg

# Read data from matrix generated with Maxwel Fortran (data from CSR format)
a = np.genfromtxt('DATA_AN_JA.data')
b = np.genfromtxt('DATA_IA.data',dtype=int)
c = np.genfromtxt('DATA_IND_AD.data')
d = np.genfromtxt('total_field_results.csv', delimiter=',', skip_header=True)

AN = a[:,0] + 1j*a[:,1] # Non diagonal values of global matrix A
AN = AN[0:-1]

JA = a[:,2].copy() # Index of the column for non diagonal terms
JA = JA.astype(dtype=int)
JA = JA - 1
JA = JA[0:-1]

IA = b.copy() # Index of starting (in row i) nondiagonal term in AN
IA = IA -1 

IND = c[:,0] + 1j*c[:,1] # Independent vector
indep_vect = IND

AD = c[:,2] + 1j*c[:,3] # Diagonal values of global matrix A

x = d[:,0].copy() # Position x
y = d[:,1].copy() # Position y


NP = len(AD)  # Number of nodes
NNZ = IA[-1]  # total number of off-diagonal entries

# Build lists to accumulate full CSR data
data = []
indices = []
indptr = [0]

# Loop over rows to insert diagonal + off-diagonals
for i in range(NP):
    row_entries = [] # Notice that row entries is initialized at each iteration

    # Diagonal entry
    row_entries.append((i, AD[i]))  # (col, value)

    # Off-diagonals: from IA[i] .. IA[i+1]-1
    for k in range(IA[i], IA[i+1]):
        row_entries.append((JA[k], AN[k]))

    # Sort by column index (important for SciPy CSR consistency)
    row_entries.sort(key=lambda x: x[0])

    # Append to global arrays
    for col, val in row_entries:
        indices.append(col)
        data.append(val)

    indptr.append(len(data))

# Convert to numpy arrays
data = np.array(data, dtype=complex)    # keep complex if needed
indices = np.array(indices, dtype=int)
indptr = np.array(indptr, dtype=int)

# Build CSR matrix
A = csr_matrix((data, indices, indptr), shape=(NP, NP))

# Solve Ax = b
H, info = bicg(A, indep_vect)

if info == 0:
    print("Converged solution:", H)
else:
    print("BiCG did not converge, info =", info)

with open('checker.csv', 'w') as f:
    f.write('X,Y,Htot_real_z,Htot_imag_z\n')
    for xi,yi,reHi,imHi in zip(x,y,np.real(H),np.imag(H)):
        f.write(f"{xi},{yi},{reHi},{imHi}\n")