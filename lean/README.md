# CutValidity: Formal Verification of Proposition 2 in Lean 4

Formal verification in Lean 4 of the cut validity results from:
> **"A Logic-Based Benders Decomposition Approach for General Lot-Sizing and Scheduling Problem with Sequence-Dependent Setup Costs"**  
> *Eyüp Ensar Işık, Z. Caner Taşkın, Semra Ağralı (2026)*

## Overview of Proposition 2

**Proposition 2.** Let $F_n$ be the set of active production decisions at iteration $n$, and let $\hat{\theta}_n = \nu(\hat{x})$ be the corresponding optimal subproblem setup cost. Under a triangular setup cost structure, the lifted cut (27) that has $LLB_n$ as the lower bound $L$ is a valid cutting plane for the Master Problem:
$$\theta \ge (\hat{\theta}_n - LLB_n) \left( \sum_{p,t \in F_n} x_{pt} \right) - (\hat{\theta}_n - LLB_n)(|F_n| - 1) + LLB_n$$

For any future master problem solution $F^*$ with preserved components $\kappa = |F_n \cap F^*|$ and omitted components $k = |F_n| - \kappa$:
$$\text{RHS} = \hat{\theta}_n - k(\hat{\theta}_n - LLB_n)$$

### Proof Structure & Verification

The paper establishes validity in four steps:
1. **Step 1 (Lower bound on $\nu$ via outgoing costs):** For any set $F$, $\nu(F) \ge \psi(F)$, where $\psi(F) = \sum_{t \in T} \left[ \sum_{j \in F_t} msc_{p_j} - \max_{j \in F_t} msc_{p_j} \right]$.
2. **Step 2 (Monotonicity under triangle inequality):** For $G = F_n \cap F^* \subseteq F^*$, metric shortcutting ensures $\nu(F^*) \ge \nu(G)$.
3. **Step 3 (Bounding decrease under job removal):** For $R = F_n \setminus G$ with $|R| = k$, $\psi(G) = \psi(F_n \setminus R) \ge \psi(F_n) - k(\psi(F_n) - LLB_n)$.
   - *Note on paper's proof:* The paper uses the heuristic intuition that "fewer jobs means a weakly smaller sum-minus-max". While $A \le B$ does not generally imply $\Delta A \le \Delta B$, our Lean formalization proves that the marginal contribution function is strictly supermodular:
     $$\Delta(S, x) = \min(c(x), \max(S \setminus \{x\}))$$
     Since $S' \subseteq S \implies \max(S' \setminus \{x\}) \le \max(S \setminus \{x\})$, the marginal gain is monotonically non-decreasing, which rigorously justifies the induction step.
4. **Step 4 (Algebraic substitution of $\hat{\theta}_n$):** Since $k \ge 1$ and $\hat{\theta}_n \ge \psi(F_n)$, $(k - 1)(\hat{\theta}_n - \psi(F_n)) \ge 0$, which yields:
   $$\psi(F_n) - k(\psi(F_n) - LLB_n) \ge \hat{\theta}_n - k(\hat{\theta}_n - LLB_n) = \text{RHS}$$

## Project Modules

- `CutValidity/PathCost.lean`:
  - `pathCost d l`: Setup cost of sequence `l`.
  - `pathCost_shortcut`: Metric shortcutting lemma (removing an intermediate job never increases setup cost under triangle inequality).
  - `pathCost_ge_sum_sub_max`: Formal verification of within-period bound $\nu(F_t) \ge \sum msc - \max msc$.

- `CutValidity/Supermodular.lean`:
  - `maxZero c S`: Maximum over finset with empty-set default 0.
  - `intraCost c S`: Baseline surrogate $\sum_{x \in S} c(x) - \max_{x \in S} c(x)$.
  - `intraCost_marginal`: Exact marginal formula $\Delta(S, x) = \min(c(x), \max(S \setminus \{x\}))$.
  - `intraCost_marginal_mono`: Formal proof of supermodularity.

- `CutValidity/Basic.lean`:
  - `cutRHS`: Definition of the affine RHS of Cut (27).
  - `cutRHS_eq_sub_m`: Identity $\text{RHS} = \hat{\theta}_n - k(\hat{\theta}_n - L)$.
  - `step4_bound`: Formal algebraic proof of Step 4.
  - `step3_induction`: Formal inductive proof of Step 3 using supermodularity.
  - `proposition2_cut_validity`: Master theorem proving unconditional validity of Cut (27) for any master solution $F^*$.
  - `cut_validity_with_GLB`: Formal verification of Proposition 1.

## Build and Verify

```bash
lake build
```
Builds cleanly with zero errors and zero `sorry`s.
