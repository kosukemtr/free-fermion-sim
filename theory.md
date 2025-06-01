# Theory of Efficient Sampling from Free Fermionic Dynamics

## Abstract

We present the theoretical background for a simulator that samples from
fermionic Gaussian states obtained by evolving Fock states under free
fermion (quadratic) Hamiltonians. Using the Majorana formalism we show
that the sampling and probability evaluation can be performed with
polynomial complexity in the number of fermionic modes.

## 1. Introduction

Free fermionic systems, characterised by quadratic Hamiltonians in
creation and annihilation operators, appear in condensed matter
physics and quantum information. These systems are exactly solvable and
their dynamics map Gaussian states to Gaussian states. A Gaussian state
is fully specified by its covariance matrix, allowing classical
simulation with resources scaling polynomially with the number of
modes. The goal of the simulator is twofold:

1. **Sample generation** – obtain bitstrings corresponding to number
   states drawn from an evolved Gaussian state.
2. **Probability evaluation** – compute the probability of a specific
   number state conditioned on the Gaussian state's covariance.

## 2. Gaussian states in the Majorana picture

For a set of $N$ fermionic modes with annihilation operators
$a_j$, we define $2N$ Hermitian Majorana operators

$$
\gamma_{2j}   = a_j + a_j^\dagger,
\qquad
\gamma_{2j+1} = -i (a_j - a_j^\dagger).
$$

A fermionic Gaussian state is completely characterised by its real,
antisymmetric covariance matrix

$$
\Gamma_{kl} = \frac{i}{2} \, \mathrm{Tr}\bigl( \rho [\gamma_k,\gamma_l] \bigr),
$$

where $\rho$ is the density matrix. Pure Gaussian states satisfy
$\Gamma^2 = -\mathbb{1}$ and so the eigenvalues of $\Gamma$ occur in
pairs $\pm i$.

A Fock basis state with occupied set $\mathcal{O}$ has block-diagonal
covariance

$$
\Gamma = \bigoplus_{j=0}^{N-1}
\begin{pmatrix}
0 & s_j \\
-s_j & 0
\end{pmatrix},
\qquad s_j = \begin{cases}
 +1 & j\in\mathcal{O},\\
 -1 & j\notin\mathcal{O}.
\end{cases}
$$

## 3. Free evolution and Bogoliubov transformations

A quadratic Hamiltonian can be written as $H = \tfrac{i}{4}
\sum_{jk} A_{jk} \gamma_j \gamma_k$ where $A$ is real and
antisymmetric. Time evolution by $H$ induces the orthogonal
transformation

$$
\Gamma \mapsto O \Gamma O^T,
\qquad O = e^{A t} \in \mathrm{SO}(2N).
$$

Likewise, a single-particle unitary $U$ acting on the mode operators
induces an orthogonal matrix in the Majorana basis via

$$
O_U = \begin{pmatrix}
\mathrm{Re}\,U & -\mathrm{Im}\,U \\
\mathrm{Im}\,U & \mathrm{Re}\,U
\end{pmatrix}.
$$

Therefore, evolving a Fock state requires only matrix multiplications
on $\Gamma$.

## 4. Sampling and probabilities

For pure Gaussian states, the probability of measuring a bitstring
$\mathbf{n} = (n_0,\ldots,n_{N-1})$ is given by a Pfaffian of a submatrix
of $\Gamma$. Concretely, let $M$ be the set of Majorana indices
corresponding to occupied modes in $\mathbf{n}$. The probability reads

$$
P(\mathbf{n}) = \mathrm{Pf}(\Gamma_M),
$$

where $\Gamma_M$ is obtained by selecting rows and columns indexed by
$M$. Computing a Pfaffian costs $\mathcal{O}(N^3)$ time using Gaussian
elimination, so enumerating all $2^N$ probabilities is exponential in
$N$, but evaluating a single outcome is polynomial.

To obtain a sample, one may either compute the full distribution when
$N$ is small or perform sequential sampling mode by mode using conditional
covariances. Both approaches only require matrix operations scaling
polynomially with $N$.

## 5. Complexity analysis

The heavy operations are matrix exponentiation to obtain $O$ (typically
$\mathcal{O}(N^3)$ using standard algorithms) and Pfaffian evaluation
(also $\mathcal{O}(N^3)$). The simulator therefore performs all tasks
— state evolution, sampling, and single-outcome probability evaluation
— with polynomial complexity in the number of modes.

## 6. Summary

Using the Majorana formalism reduces the dynamics of free fermions to
real orthogonal transformations of covariance matrices. This approach
allows efficient classical simulation of evolved Gaussian states.
Sampling and computing probabilities of number states require only
linear algebra and scale polynomially with the system size.
