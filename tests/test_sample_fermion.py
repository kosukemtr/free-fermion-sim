import unittest
import numpy as np

from fermion_sampler.sample_fermion import (
    majorana_hamiltonian,
    orthogonal_from_hamiltonian,
    covariance_from_occupation,
    evolve_covariance,
    orbital_to_majorana_rotation,
    sample_fermion_state,
    get_bitstring_probs_from_covariance,
    sample_evolved_state,
)

class TestFermionSampler(unittest.TestCase):
    def test_majorana_hamiltonian_antisymmetric(self):
        rng = np.random.default_rng(123)
        x = rng.normal(size=(4, 4))
        A = majorana_hamiltonian(x)
        self.assertTrue(np.allclose(A, -A.T))

    def test_orthogonal_from_hamiltonian(self):
        rng = np.random.default_rng(456)
        x = rng.normal(size=(6, 6))
        A = majorana_hamiltonian(x)
        O = orthogonal_from_hamiltonian(A)
        eye = np.eye(O.shape[0])
        self.assertTrue(np.allclose(O @ O.T, eye, atol=1e-12))

    def test_orbital_to_majorana_rotation_orthogonal(self):
        rng = np.random.default_rng(789)
        X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        U, _ = np.linalg.qr(X)
        O = orbital_to_majorana_rotation(U)
        eye = np.eye(O.shape[0])
        self.assertTrue(np.allclose(O @ O.T, eye, atol=1e-12))

    def test_evolve_covariance_preserves_antisymmetry(self):
        rng = np.random.default_rng(101)
        x = rng.normal(size=(4, 4))
        A = majorana_hamiltonian(x)
        O = orthogonal_from_hamiltonian(A)
        Gamma = A  # any antisymmetric matrix works as a covariance analogue
        evolved = evolve_covariance(Gamma, O)
        self.assertTrue(np.allclose(evolved, -evolved.T, atol=1e-12))

    def test_covariance_from_occupation_blocks(self):
        N = 3
        occ = [0, 2]
        Gamma = covariance_from_occupation(occ, N)
        for j in range(N):
            block = Gamma[2 * j : 2 * j + 2, 2 * j : 2 * j + 2]
            expected = 1 if j in occ else -1
            ref = np.array([[0, expected], [-expected, 0]])
            self.assertTrue(np.array_equal(block, ref))

    def test_sampling_matches_probabilities(self):
        N = 2
        occ = [1]
        Gamma = covariance_from_occupation(occ, N)
        probs = get_bitstring_probs_from_covariance(Gamma)
        counts = {i: 0 for i in range(len(probs))}
        for _ in range(2000):
            bitstr = sample_fermion_state(Gamma)
            idx = int(bitstr, 2)
            counts[idx] += 1
        freq = np.array([counts[i] / 2000 for i in range(len(probs))])
        self.assertTrue(np.allclose(freq, probs, atol=0.1))

    def test_sample_evolved_state_length(self):
        N = 4
        occ = [0, 1]
        U = np.eye(N)
        out = sample_evolved_state(occ, U)
        self.assertIsInstance(out, str)
        self.assertEqual(len(out), N)


if __name__ == "__main__":
    unittest.main()
