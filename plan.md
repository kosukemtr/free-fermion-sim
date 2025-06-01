# Implementation Plan

This document outlines the steps required to realise a simulator for
free fermion dynamics capable of sampling number states and evaluating
their probabilities.

## 1. Core Library

- **Majorana utilities** – functions converting between orbital and
  Majorana representations and constructing antisymmetric Hamiltonians.
- **Covariance manipulation** – routines to generate covariance matrices
  for Fock states, apply orthogonal transformations and compute
  probabilities via Pfaffians.
- **Sampling** – sequential or direct sampling algorithms operating on
  covariance matrices.

## 2. API Design

Create a module `fermion_sampler` exposing:

- `majorana_hamiltonian(x)` – build antisymmetric generator from any
  real matrix.
- `orthogonal_from_hamiltonian(A)` – obtain `exp(A)`.
- `covariance_from_occupation(occ, N)` – covariance of a Fock state.
- `evolve_covariance(Gamma, O)` – apply orthogonal evolution.
- `orbital_to_majorana_rotation(U)` – convert a unitary to a real
  orthogonal matrix.
- `get_bitstring_probs_from_covariance(Gamma)` – compute the probability
  vector.
- `sample_fermion_state(Gamma)` – draw a single bitstring sample.
- `sample_evolved_state(occ, U)` – convenience routine performing all
  steps from initial occupation to sampling after evolution.

## 3. Testing

- Unit tests validating antisymmetry, orthogonality and sampling
  behaviour.
- Numerical checks comparing empirical frequencies against theoretical
  probabilities for small systems.

## 4. Extensions

- Support mixed Gaussian states using covariance and displacement.
- Implement more efficient sampling algorithms for large mode numbers
  (e.g. conditional sampling).
- Interface with external quantum frameworks (Cirq, QURI Parts) to
  construct circuits corresponding to the simulated dynamics.

## 5. Example Usage

Provide notebooks and scripts demonstrating evolution of simple Fock
states, probability evaluation and sampling.
