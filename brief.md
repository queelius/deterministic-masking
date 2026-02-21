# Information Recovery under Deterministic Masking in Exponential Series Systems

## Status: Draft (proof complete, simulation validated)

## Thesis

For exponential series systems, a deterministic masking mechanism that
injectively maps each failure cause to a unique candidate set recovers the
*complete-identification* Fisher information matrix, regardless of masking
cardinality. Information loss from masking is a property of the masking
*mechanism*, not the masking *cardinality*.

## Core Result

**Proposition (Information Recovery under Deterministic Masking).**
For an m-component exponential series system with 1 <= w <= m-1, suppose
the masking mechanism is *separating-deterministic*: there exists an
injection phi: {1,...,m} -> {C subset {1,...,m} : |C| = w} satisfying
k in phi(k) and P(C = phi(k) | K = k) = 1 for each k. Then the Fisher
information matrix equals that of complete identification (w = 1):

    I(lambda | phi)_{jk} = delta_{jk} / (Lambda * lambda_j)

In particular, for w = m-1 the candidate sets {1,...,m}\{j} are in
bijection with the m components, and any derangement sigma defines such
a mechanism via phi(k) = {1,...,m}\{sigma(k)}. Derangements exist for
all m >= 2.

## Proof Sketch

1. Injectivity of phi means phi^{-1}(C) recovers K with certainty
2. The joint density marginalizes to f(t, c | lambda) = lambda_{phi^{-1}(c)} * exp(-Lambda*t)
3. This has the same functional form as the unmasked density f(t, k | lambda) = lambda_k * exp(-Lambda*t)
4. Same functional form => same second derivatives => same FIM
5. Key ingredient: K independent of S (exponential-specific, Prop 4.3 in parent paper)

## Key Implications

1. **Mechanism vs cardinality**: The w=m-1 uniform masking FIM is *not*
   a worst-case bound over all masking mechanisms at that cardinality.
   Deterministic masking at w=m-1 matches w=1 exactly.

2. **I(K;C) vs FIM(lambda)**: Diagnostic mutual information (how well C
   identifies the failed component) and Fisher information (how precisely
   parameters can be estimated) are *distinct* quantities. Uniform masking
   minimizes I(K;C), but deterministic masking has I(K;C) = H(K) (maximum)
   and also achieves maximum FIM.

3. **C2 violation**: The separating-deterministic mechanism violates the
   symmetry condition C2, since P(C = phi(k) | K = k) = 1 but
   P(C = phi(k) | K = k') = 0 for k' != k with k' in phi(k). The C2
   assumption is what costs information.

4. **Diagnostic design**: A well-designed diagnostic that consistently maps
   each failure mode to a unique signature incurs *no information loss
   whatsoever*. This has practical implications for test design.

## Numerical Verification

Monte Carlo simulation (R = 100,000) comparing three mechanisms on
identical (S, K) data for m = 3, 5:

| Mechanism          | Var(lambda_1) * n (m=3, lambda=(1,3,5)) |
|--------------------|----------------------------------------|
| w=1 (no masking)   | 9.00 (theory), 9.02 (empirical)        |
| Deterministic w=2  | 9.00 (theory), 9.02 (empirical)        |
| Uniform w=2        | 81.00 (theory), 80.79 (empirical)      |

Empirical difference between w=1 and deterministic: **0.0000%** across
all sample sizes (n = 50 to 2000) and both configurations.

## Available Assets

From parent paper repo (expo-masked-fim):
- `research/validate_deterministic_masking.R` — full simulation script
- `paper/fig_deterministic_masking.pdf` — two-panel comparison figure
- `research/verify_masking_optimality.R` — verifies uniform is NOT optimal

## Extensions

1. **Perturbation analysis**: How does the FIM degrade continuously as the
   masking mechanism moves from deterministic toward uniform? What is the
   functional form of the information loss Delta(epsilon) for a mechanism
   that is epsilon-close to deterministic?

2. **Non-exponential distributions**: Does the result extend to Weibull or
   other distributions where K is NOT independent of S? Conjecture: the
   result holds approximately when hazard rates vary slowly.

3. **Optimal mechanism design**: Among all masking mechanisms at cardinality w,
   which one maximizes det(I)? Is the answer always deterministic (when
   m <= C(m,w), so injection exists)?

4. **Practical diagnostic design**: When can real diagnostic systems achieve
   near-deterministic masking? Implications for test procedure design.

## Connection to Parent Paper

This result was originally Proposition 4.7 in expo-masked-fim but extracted
to keep that paper at 33 pages. The parent paper retains the conceptual
point (I(K;C) != FIM(lambda)) as a remark following Proposition 2.6,
without the formal proposition or proof.

## Target Venue

- **Primary**: Statistics & Probability Letters (short communication, 6 page limit)
- **Alternative**: IEEE Transactions on Reliability (short paper)
- **Alternative**: The American Statistician (short note, educational focus)

## Dependencies

- K independent of S for exponential distributions (Prop 4.3 in parent paper)
- FIM formula for general masking mechanisms
- Parent paper: expo-masked-fim (Towell 2025)
- Masters thesis: towell2023reliability
