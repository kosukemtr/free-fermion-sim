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

For a set of \(N\) fermionic modes with annihilation operators
\(a_j\), we define \(2N\) Hermitian Majorana operators

\[
\gamma_{2j}   \;=\; a_j + a_j^\dagger,
\qquad
\gamma_{2j+1} \;=\; -i \bigl(a_j - a_j^\dagger\bigr).
\]

A fermionic Gaussian state is completely characterised by its real,
antisymmetric covariance matrix

\[
\Gamma_{kl} \;=\; \frac{i}{2}\,
\operatorname{Tr}\!\bigl(\rho\,[\gamma_k,\gamma_l]\bigr),
\]

where \(\rho\) is the density matrix.  Pure Gaussian states satisfy
\(\Gamma^2 = -\mathbb{1}\), so the eigenvalues of \(\Gamma\) occur in
pairs \(\pm i\).

A Fock basis state with occupied set \(\mathcal{O}\) has block-diagonal
covariance

\[
\Gamma \;=\; \bigoplus_{j=0}^{N-1}
\begin{pmatrix}
0 & s_j \\
-s_j & 0
\end{pmatrix},
\qquad
s_j \;=\;
\begin{cases}
 +1 & j\in\mathcal{O},\\[4pt]
 -1 & j\notin\mathcal{O}.
\end{cases}
\]

## 3. Free evolution and Bogoliubov transformations

A quadratic Hamiltonian can be written as
\[
H \;=\; \tfrac{i}{4}\sum_{jk} A_{jk}\,\gamma_j\gamma_k,
\]
with \(A\) real and antisymmetric.  
Time evolution by \(H\) induces the orthogonal transformation

\[
\Gamma \;\longmapsto\; O\,\Gamma\,O^{\mathsf T},
\qquad
O \;=\; e^{A t} \in \mathrm{SO}(2N).
\]

Likewise, a single-particle unitary \(U\) acting on the mode operators
induces an orthogonal matrix in the Majorana basis via

\[
O_U \;=\;
\begin{pmatrix}
\operatorname{Re}U & -\operatorname{Im}U \\
\operatorname{Im}U & \phantom{-}\operatorname{Re}U
\end{pmatrix}.
\]

Evolving a Fock state therefore requires only matrix multiplications on
\(\Gamma\).

## 4. Sampling and probabilities

For pure Gaussian states, the probability of measuring a bitstring
\(\mathbf n = (n_0,\ldots,n_{N-1})\) is given by a Pfaffian of a
sub-matrix of \(\Gamma\).  Let \(M\) be the set of Majorana indices
corresponding to occupied modes in \(\mathbf n\).  Then

\[
P(\mathbf n) \;=\; \operatorname{Pf}\bigl(\Gamma_M\bigr),
\]

where \(\Gamma_M\) is obtained by selecting rows and columns indexed by
\(M\).  Because \(\Gamma\) is a valid covariance matrix, all such
Pfaffians are real numbers strictly between \(0\) and \(1\), and

\[
\sum_M \operatorname{Pf}(\Gamma_M) \;=\; 1.
\]

Computing a Pfaffian costs \(\mathcal{O}(N^3)\) using Gaussian
elimination, so enumerating all \(2^N\) probabilities is exponential in
\(N\).  Evaluating a single outcome, however, is polynomial.

To obtain a sample one may either compute the full distribution when
\(N\) is small or perform sequential sampling mode by mode using
conditional covariances.  Both approaches rely only on matrix
operations scaling polynomially with \(N\).

### 4.1 Sequential sampling

Sequential sampling constructs a bitstring one mode at a time using
conditional probabilities derived from the covariance matrix.  Let
\(\Gamma^{(j)}\) denote the covariance conditioned on the outcomes for
modes \(0,\ldots,j\!-\!1\).  The probability that mode \(j\) is occupied
is

\[
p_j \;=\; \frac{1+\Gamma^{(j)}_{2j,\,2j+1}}{2}.
\]

A Bernoulli random variable \(n_j\in\{0,1\}\) is drawn with that
probability.  To update the state we partition the conditioned
covariance matrix as

\[
\Gamma^{(j)} \;=\;
\begin{pmatrix}
  A_j & B_j \\
 -B_j^{\mathsf T} & C_j
\end{pmatrix},
\]

where \(A_j\) is the \(2\times2\) block for mode \(j\).
Denote by \(J_{n_j}\) the covariance of a Fock state with occupation
\(n_j\),

\[
J_{n_j} \;=\;
\begin{pmatrix}
 0 & s_j \\
 -s_j & 0
\end{pmatrix},
\qquad
s_j \;=\; (-1)^{1-n_j}.
\]

The conditional covariance for the remaining modes is obtained via the
Schur complement

\[
\Gamma^{(j+1)} \;=\;
C_j \;-\;
B_j^{\mathsf T}\bigl(A_j - J_{n_j}\bigr)^{-1} B_j.
\]

Iterating this procedure for all modes yields the bitstring
\((n_0,\ldots,n_{N-1})\).  Multiplying the conditional probabilities at
each step recovers the joint probability, which equals
\(\operatorname{Pf}(\Gamma_M)\) with \(M\) the set of occupied modes.
The overall cost is \(\mathcal{O}(N^3)\), dominated by solving the
linear systems in the Schur complements.

## 5. Complexity analysis

The dominant operations are matrix exponentiation to obtain \(O\)
(generally \(\mathcal{O}(N^3)\)) and Pfaffian evaluation
(\(\mathcal{O}(N^3)\)).  Consequently, state evolution, sampling, and
single-outcome probability evaluation all run in polynomial time.

## 6. Summary

Using the Majorana formalism reduces the dynamics of free fermions to
real orthogonal transformations of covariance matrices, enabling
efficient classical simulation of evolved Gaussian states.  Sampling
and computing probabilities of number states require only linear
algebra and scale polynomially with system size.
