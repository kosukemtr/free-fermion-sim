import numpy as np
from scipy.linalg import expm

# ----------------------------------------------------------------------
#  Core helpers
# ----------------------------------------------------------------------

def majorana_hamiltonian(x: np.ndarray) -> np.ndarray:
    """Return the antisymmetric Majorana matrix  A = x − xᵀ."""
    return x - x.T                                   # ensure antisymmetry


def orthogonal_from_hamiltonian(A: np.ndarray) -> np.ndarray:
    """Compute the real orthogonal matrix  O = exp(−A)  for a 2N×2N antisymmetric A."""
    return expm(-A)


def covariance_from_occupation(occ, N: int) -> np.ndarray:
    """
    Build the 2N×2N Majorana covariance matrix Γ corresponding to
    the Fock state with occupied mode indices in *occ*.
    """
    Γ = np.zeros((2 * N, 2 * N))
    for p in occ:
        Γ[2 * p,     2 * p + 1] =  1
        Γ[2 * p + 1, 2 * p    ] = -1
    return Γ


def evolve_covariance(Γ: np.ndarray, O: np.ndarray) -> np.ndarray:
    """Time-evolve a covariance matrix:  Γ ↦ O Γ Oᵀ ."""
    return O @ Γ @ O.T


def orbital_to_majorana_rotation(U: np.ndarray) -> np.ndarray:
    """
    Convert an N×N **unitary** orbital-rotation matrix U into the 2N×2N real
    orthogonal matrix O acting on Majorana operators, so that γ′ = O γ.
    """
    if U.ndim != 2 or U.shape[0] != U.shape[1]:
        raise ValueError("U must be a square matrix")

    O = np.block([[ U.real, -U.imag],
                  [ U.imag,  U.real]])

    # sanity: O should be orthogonal if U is unitary
    if not np.allclose(O @ O.T, np.eye(2 * U.shape[0]), atol=1e-12):
        raise ValueError("Result is not orthogonal—input U may not be unitary")
    return O


# ----------------------------------------------------------------------
#  Sampling and probability utilities
# ----------------------------------------------------------------------

def sample_fermion_state(Γ: np.ndarray) -> str:
    """
    Draw a single exact sample from the Gaussian (Slater-determinant) state
    specified by its Majorana covariance matrix Γ, and return the result as
    a *bit-string*.

    The returned string is **big-endian**:
        • The left-most character is mode N−1  
        • The right-most character is mode 0

    Example for N = 4:  '1001'  represents |1 0 0 1⟩ with modes (3 2 1 0).
    """
    N = Γ.shape[0] // 2
    # Build particle-sector correlator  C = (I + iJ Γ)/2
    J = np.kron(np.eye(N), np.array([[0, -1], [1, 0]]))
    C = 0.5 * (np.eye(2 * N) + 1j * J @ Γ)
    C = C[:N, :N].copy()                   # restrict to particle sector

    bits_le = []                           # collect bits little-endian first
    rng = np.random.default_rng()
    for _ in range(N):
        p_occ = rng.random() < C[0, 0].real
        bits_le.append(int(p_occ))

        if p_occ:
            # project out the occupied mode (Sherman-Morrison update)
            v = C[:, 0].copy()
            C -= np.outer(v, v.conj()) / C[0, 0]

        # Remove the measured mode
        C = C[1:, 1:]

    # Convert to big-endian string
    return ''.join(str(b) for b in reversed(bits_le))


def sample_evolved_state(occ, U: np.ndarray) -> str:
    """
    Take an initial occupation list *occ*, evolve with orbital unitary U,
    and return one measurement outcome as a bit-string (big-endian).
    """
    N = len(occ)
    O = orbital_to_majorana_rotation(U)
    Γ = covariance_from_occupation(occ, N)
    Γ_evolved = evolve_covariance(Γ, O)
    return sample_fermion_state(Γ_evolved)


def get_bitstring_probs_from_covariance(Γ: np.ndarray) -> np.ndarray:
    """
    Enumerate the full probability vector P of length 2**N for measuring
    every occupation-basis bit-string, given a covariance matrix Γ.

    *Intended for debugging / validation of small systems (N ≲ 12–14).*
    """
    N = Γ.shape[0] // 2
    if Γ.shape != (2 * N, 2 * N):
        raise ValueError("Γ must be 2N×2N")

    # Particle-sector correlator
    J = np.kron(np.eye(N), np.array([[0, -1], [1, 0]]))
    C_root = 0.5 * (np.eye(2 * N) + 1j * J @ Γ)
    C_root = C_root[:N, :N]

    def recurse(C_sub):
        m = C_sub.shape[0]
        if m == 0:
            return {(): 1.0}

        p00 = C_sub[0, 0].real

        # branch: unoccupied
        dist0 = recurse(C_sub[1:, 1:])

        # branch: occupied
        v = C_sub[:, 0].copy()
        C1 = C_sub - np.outer(v, v.conj()) / C_sub[0, 0]
        dist1 = recurse(C1[1:, 1:])

        out = {}
        for bits, p in dist0.items():
            out[(0,) + bits] = (1 - p00) * p
        for bits, p in dist1.items():
            out[(1,) + bits] = p00 * p
        return out

    dist = recurse(C_root)

    probs = np.zeros(1 << N, dtype=float)
    for bits, p in dist.items():
        idx = sum((b << i) for i, b in enumerate(bits))  # little-endian → int
        probs[idx] = p.real

    if not np.allclose(probs.sum(), 1.0, atol=1e-12):
        raise RuntimeError("Probabilities do not sum to 1")
    return probs