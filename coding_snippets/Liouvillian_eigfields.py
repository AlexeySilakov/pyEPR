import numpy as np

###############################################################################
# 1. Spin matrices for S=1/2 and I=3/2
###############################################################################

def spin_matrices(S):
    """Return Sx, Sy, Sz for spin S."""
    dim = int(2*S+1)
    m = np.arange(S, -S-1, -1)  # m = S, S-1, ..., -S
    Splus = np.zeros((dim,dim))
    Sminus = np.zeros((dim,dim))

    for i, mi in enumerate(m):
        for j, mj in enumerate(m):
            if mj == mi-1:  # lowering from mi to mj
                Sminus[i,j] = np.sqrt(S*(S+1) - mj*(mj+1))
            if mj == mi+1:  # raising
                Splus[i,j] = np.sqrt(S*(S+1) - mj*(mj-1))

    Sx = 0.5*(Splus + Sminus)
    Sy = -0.5j*(Splus - Sminus)
    Sz = np.diag(m)
    return Sx, Sy, Sz

# electron (S=1/2)
Sx, Sy, Sz = spin_matrices(1/2)

# nucleus (I=3/2)
Ix, Iy, Iz = spin_matrices(3/2)

# Tensor products for total Hamiltonian
def kron(a,b): return np.kron(a,b)

# Build total spin operators for the system
Sx_tot = kron(Sx, np.eye(4))
Sy_tot = kron(Sy, np.eye(4))
Sz_tot = kron(Sz, np.eye(4))

Ix_tot = kron(np.eye(2), Ix)
Iy_tot = kron(np.eye(2), Iy)
Iz_tot = kron(np.eye(2), Iz)

###############################################################################
# 2. Hamiltonian for EPR: H = μB g B Sz + A (S · I)
###############################################################################

muB = 13.99624555  # GHz/T (approx)
g = 2.0023
A = 0.10            # 100 MHz hyperfine (example)
omega = 9.5         # microwave frequency (GHz)

def Hamiltonian(B):
    """Returns 8x8 Hamiltonian for field B (Tesla)."""
    Hzeeman = muB * g * B * Sz_tot
    Hhf = A * (Sx_tot@Ix_tot + Sy_tot@Iy_tot + Sz_tot@Iz_tot)
    return Hzeeman + Hhf

###############################################################################
# 3. Liouvillian superoperator L = H⊗I - I⊗H^T
###############################################################################

def liouvillian(H):
    dim = H.shape[0]
    return np.kron(H, np.eye(dim)) - np.kron(np.eye(dim), H.T)

###############################################################################
# 4. Solve eigenfield condition: L ρ = i ℏ ω ρ
###############################################################################

def solve_eigenfield(B):
    H = Hamiltonian(B)
    L = liouvillian(H)
    # frequency term for GHz units: ℏ = 1 here (units consistent)
    M = L - 1j*omega*np.eye(L.shape[0])
    # Solve eigenproblem
    vals, vecs = np.linalg.eig(M)
    # Return eigenvectors with very small eigenvalues
    idx = np.argsort(np.abs(vals))  # smallest magnitude = resonance
    return vals[idx], vecs[:,idx], H


###############################################################################
# 5. Identify dominant transition from Liouvillian eigenvector
###############################################################################

def identify_transition(H, rho_vec):
    """Return (m,n) transition and dominant |mS,mI> states."""
    dim = H.shape[0]

    # Unvectorize
    R = rho_vec.reshape((dim,dim))

    # Diagonalize Hamiltonian
    E, U = np.linalg.eigh(H)

    # Transform R into energy basis
    R_tilde = U.conj().T @ R @ U

    # Find largest off-diagonal element
    mask = np.ones_like(R_tilde, dtype=bool)
    np.fill_diagonal(mask, 0)
    i, j = np.unravel_index(np.argmax(np.abs(R_tilde[mask])), R_tilde.shape)

    # Return results
    return (i, j), E, R_tilde


###############################################################################
# 6. Demonstration: scan field and identify the first resonance
###############################################################################

B_test = 0.34   # Tesla, near X-band resonance
vals, vecs, H = solve_eigenfield(B_test)

# Take the eigenvector with smallest |eigenvalue|
rho_vec = vecs[:,0]

(m,n), E, R_tilde = identify_transition(H, rho_vec)

print("=== Resonance analysis at B = %.4f T ===" % B_test)
print("Dominant transition: level m =", m, ", level n =", n)
print("Energy difference: %.4f GHz" % (E[m] - E[n]))
print("Expected microwave frequency: %.4f GHz" % omega)

###############################################################################
# 7. Convert (m,n) levels into (mS, mI) quantum numbers
###############################################################################

# Basis ordering used: |mS, mI> with mS = +1/2,-1/2 and mI= +3/2,+1/2,-1/2,-3/2
mS_vals = [+0.5, +0.5, +0.5, +0.5, -0.5, -0.5, -0.5, -0.5]
mI_vals = [+1.5, +0.5, -0.5, -1.5, +1.5, +0.5, -0.5, -1.5]

print("\nTransition quantum numbers:")
print("(mS,mI) for state m =", (mS_vals[m], mI_vals[m]))
print("(mS,mI) for state n =", (mS_vals[n], mI_vals[n]))
