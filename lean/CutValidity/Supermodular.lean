import Mathlib.Data.Finset.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Data.Finset.Lattice.Fold
import Mathlib.Tactic

set_option linter.style.header false
set_option linter.style.longLine false
set_option linter.unusedSectionVars false
set_option linter.deprecated false
set_option linter.unusedDecidableInType false

open Finset

variable {α : Type*} [DecidableEq α]

/-!
# Supermodularity and Marginal Contribution of Intra-Period Setup Costs

This module proves that the baseline intra-period surrogate function
`intraCost(S) = ∑_{x ∈ S} c(x) - max_{x ∈ S} c(x)` is supermodular.
Specifically:
1. `intraCost_marginal`: The marginal contribution of job `x` to `intraCost(S)`
   equals `min(c(x), max_{y ∈ S \ {x}} c(y))`.
2. `intraCost_marginal_mono`: If `S' ⊆ S` and `x ∈ S'`, then
   `intraCost(S') - intraCost(S' \ {x}) ≤ intraCost(S) - intraCost(S \ {x})`.
3. This rigorously justifies Step 3 in the proof of Proposition 2 in
   Işık, Taşkın, Ağralı (2026), resolving the heuristic leap in the original paper.
-/

/-- Maximum of `c` over `S`, defaulting to 0 for the empty set. -/
noncomputable def maxZero (c : α → ℝ) (S : Finset α) : ℝ :=
  if h : S.Nonempty then S.sup' h c else 0

@[simp] lemma maxZero_empty (c : α → ℝ) : maxZero c ∅ = 0 := by
  simp [maxZero]

lemma maxZero_nonneg (c : α → ℝ) (hc : ∀ a, 0 ≤ c a) (S : Finset α) :
    0 ≤ maxZero c S := by
  by_cases h : S.Nonempty
  · dsimp [maxZero]
    rw [dif_pos h]
    obtain ⟨a, ha⟩ := h
    have : c a ≤ S.sup' ⟨a, ha⟩ c := le_sup' c ha
    have : 0 ≤ c a := hc a
    linarith
  · dsimp [maxZero]
    rw [dif_neg h]

lemma maxZero_le_of_subset (c : α → ℝ) (hc : ∀ a, 0 ≤ c a) {A B : Finset α} (hAB : A ⊆ B) :
    maxZero c A ≤ maxZero c B := by
  by_cases hA : A.Nonempty
  · have hB : B.Nonempty := Nonempty.mono hAB hA
    dsimp [maxZero]
    rw [dif_pos hA, dif_pos hB]
    rw [sup'_le_iff]
    intro b hb
    exact le_sup' c (hAB hb)
  · dsimp [maxZero]
    rw [dif_neg hA]
    exact maxZero_nonneg c hc B

lemma maxZero_insert (c : α → ℝ) (hc : ∀ a, 0 ≤ c a) (x : α) (S : Finset α) :
    maxZero c (insert x S) = max (c x) (maxZero c S) := by
  have h_ins : (insert x S).Nonempty := insert_nonempty x S
  by_cases hS : S.Nonempty
  · dsimp [maxZero]
    rw [dif_pos h_ins, dif_pos hS]
    rw [sup'_insert (s := S) (b := x) (f := c)]
  · have hS_eq : S = ∅ := not_nonempty_iff_eq_empty.mp hS
    subst hS_eq
    dsimp [maxZero]
    exact (max_eq_left (hc x)).symm

/-- Algebraic identity: `max a b - b = max (a - b) 0`. -/
lemma max_sub_right (a b : ℝ) : max a b - b = max (a - b) 0 := by
  rcases le_total a b with h | h
  · rw [max_eq_right h, max_eq_right (by linarith)]
    ring
  · rw [max_eq_left h, max_eq_left (by linarith)]

/-- Algebraic identity: `a - max (a - b) 0 = min a b`. -/
lemma sub_max_sub_zero (a b : ℝ) : a - max (a - b) 0 = min a b := by
  rcases le_total a b with h | h
  · rw [max_eq_right (by linarith), sub_zero, min_eq_left h]
  · rw [max_eq_left (by linarith), min_eq_right h]
    ring

/-- The within-period lower bound function: sum of costs minus maximum. -/
noncomputable def intraCost (c : α → ℝ) (S : Finset α) : ℝ :=
  S.sum c - maxZero c S

/-- The exact marginal contribution of an element `x` to `intraCost S`. -/
theorem intraCost_marginal (c : α → ℝ) (hc : ∀ a, 0 ≤ c a) (x : α) (S : Finset α) (hx : x ∈ S) :
    intraCost c S - intraCost c (S \ {x}) = min (c x) (maxZero c (S \ {x})) := by
  have h_ins : insert x (S \ {x}) = S := insert_sdiff_self_of_mem hx
  have h_not : x ∉ S \ {x} := by simp
  have h_sum : S.sum c = (S \ {x}).sum c + c x := by
    conv_lhs => rw [← h_ins]
    rw [sum_insert h_not, add_comm]
  dsimp [intraCost]
  rw [h_sum]
  have h_max : maxZero c S = max (c x) (maxZero c (S \ {x})) := by
    conv_lhs => rw [← h_ins]
    exact maxZero_insert c hc x (S \ {x})
  rw [h_max]
  have : ((S \ {x}).sum c + c x - max (c x) (maxZero c (S \ {x}))) -
         ((S \ {x}).sum c - maxZero c (S \ {x})) =
         c x - (max (c x) (maxZero c (S \ {x})) - maxZero c (S \ {x})) := by ring
  rw [this, max_sub_right (c x) (maxZero c (S \ {x})), sub_max_sub_zero (c x) (maxZero c (S \ {x}))]

/-- Monotonicity of marginal contribution (supermodularity):
    Removing an element from a smaller set reduces `intraCost` by AT MOST
    removing it from a larger set. -/
theorem intraCost_marginal_mono (c : α → ℝ) (hc : ∀ a, 0 ≤ c a)
    {S' S : Finset α} (h_sub : S' ⊆ S) (x : α) (hx : x ∈ S') :
    intraCost c S' - intraCost c (S' \ {x}) ≤ intraCost c S - intraCost c (S \ {x}) := by
  have hxS : x ∈ S := h_sub hx
  rw [intraCost_marginal c hc x S' hx, intraCost_marginal c hc x S hxS]
  have h_sub_diff : S' \ {x} ⊆ S \ {x} := sdiff_subset_sdiff h_sub (subset_refl _)
  have h_max_le := maxZero_le_of_subset c hc h_sub_diff
  exact min_le_min le_rfl h_max_le
