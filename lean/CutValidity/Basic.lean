import Mathlib.Data.Finset.Basic
import Mathlib.Data.Finset.Card
import Mathlib.Tactic
import CutValidity.Supermodular

set_option linter.style.header false
set_option linter.style.longLine false
set_option linter.unusedSectionVars false
set_option linter.deprecated false

open Finset

variable {α : Type*} [DecidableEq α]

/-!
# Formal Verification of Lifted Logic-Based Benders Cuts

This module formally verifies the validity of the lifted Benders optimality cut (27)
from the paper:
"A Logic-Based Benders Decomposition Approach for General Lot-Sizing and Scheduling Problem
with Sequence-Dependent Setup Costs" (Işık, Taşkın, Ağralı, 2026).

We verify:
1. Proposition 1: Cut validity under a Global Lower Bound (GLB).
2. Proposition 2: Cut validity under the Local Lower Bound (LLB_n), fully formalizing
   all 4 steps of the paper's proof:
   - Step 1: Subproblem lower bounding via outgoing costs (ψ ≤ ν).
   - Step 2: Monotonicity under set inclusion / metric shortcutting (ν(G) ≤ ν(F*)).
   - Step 3: Induction on the number of removed jobs m = |R| using supermodularity.
   - Step 4: Algebraic replacement of ψ(Fn) by θ_hat for m ≥ 1.
   - Master Theorem: Unconditional validity for all master problem solutions F*.
-/

/-- The Right-Hand Side of the lifted Benders optimality cut (27)
    as a function of the number of preserved components κ. -/
def cutRHS (card_Fn : ℕ) (θ_hat : ℝ) (L : ℝ) (κ : ℕ) : ℝ :=
  (θ_hat - L) * (κ : ℝ) - (θ_hat - L) * ((card_Fn : ℝ) - 1) + L

/-- Alternative simplified expression for cutRHS: (θ_hat - L) * (κ - card_Fn + 1) + L -/
lemma cutRHS_eq (card_Fn : ℕ) (θ_hat L : ℝ) (κ : ℕ) :
    cutRHS card_Fn θ_hat L κ = (θ_hat - L) * ((κ : ℝ) - (card_Fn : ℝ) + 1) + L := by
  dsimp [cutRHS]
  ring

/-- Exact representation of cutRHS in terms of the number of removed components m = card_Fn - κ. -/
lemma cutRHS_eq_sub_m (card_Fn : ℕ) (θ_hat L : ℝ) (κ : ℕ) (hκ : κ ≤ card_Fn) :
    cutRHS card_Fn θ_hat L κ = θ_hat - ((card_Fn - κ : ℕ) : ℝ) * (θ_hat - L) := by
  dsimp [cutRHS]
  have hcast : ((card_Fn - κ : ℕ) : ℝ) = (card_Fn : ℝ) - (κ : ℝ) := Nat.cast_sub hκ
  rw [hcast]
  ring

/-- When κ = card_Fn, the cut RHS evaluates exactly to θ_hat. -/
lemma cutRHS_at_card (card_Fn : ℕ) (θ_hat L : ℝ) :
    cutRHS card_Fn θ_hat L card_Fn = θ_hat := by
  rw [cutRHS_eq]
  ring

/-- When κ = card_Fn - 1 (with card_Fn ≥ 1), the cut RHS evaluates exactly to L. -/
lemma cutRHS_at_card_sub_one (card_Fn : ℕ) (hFn : 1 ≤ card_Fn) (θ_hat L : ℝ) :
    cutRHS card_Fn θ_hat L (card_Fn - 1) = L := by
  rw [cutRHS_eq]
  have hcast : ((card_Fn - 1 : ℕ) : ℝ) = (card_Fn : ℝ) - 1 := Nat.cast_sub hFn
  rw [hcast]
  ring

/-- When κ ≤ card_Fn - 2 and θ_hat ≥ L, the cut RHS is bounded above by 2 * L - θ_hat. -/
lemma cutRHS_le_two_L_sub_theta (card_Fn : ℕ) (θ_hat L : ℝ) (κ : ℕ)
    (h_pos : 0 ≤ θ_hat - L) (hκ : κ + 2 ≤ card_Fn) :
    cutRHS card_Fn θ_hat L κ ≤ 2 * L - θ_hat := by
  rw [cutRHS_eq]
  have h_le : (κ : ℝ) - (card_Fn : ℝ) + 1 ≤ -1 := by
    have h_nat : (κ : ℝ) + 2 ≤ (card_Fn : ℝ) := by exact_mod_cast hκ
    linarith
  have h_mul : (θ_hat - L) * ((κ : ℝ) - (card_Fn : ℝ) + 1) ≤ (θ_hat - L) * (-1) := by
    nlinarith
  linarith

/-- Combinatorial Lemma: If |Fn ∩ F*| = |Fn| - 1, then exactly one element j ∈ Fn
    is omitted from F*, which means Fn \ {j} ⊆ F*. -/
lemma exists_singleton_diff_of_card_inter_eq_sub_one
    (Fn F_star : Finset α) (hFn : 1 ≤ Fn.card)
    (hκ : (Fn ∩ F_star).card = Fn.card - 1) :
    ∃ j ∈ Fn, Fn \ {j} ⊆ F_star := by
  have h_sub : Fn ∩ F_star ⊆ Fn := inter_subset_left
  have h_diff_card : (Fn \ (Fn ∩ F_star)).card = 1 := by
    rw [card_sdiff_of_subset h_sub, hκ]
    omega
  obtain ⟨j, hj⟩ := card_eq_one.mp h_diff_card
  have hj_in_diff : j ∈ Fn \ (Fn ∩ F_star) := by
    rw [hj]
    exact mem_singleton_self j
  have hj_in_Fn : j ∈ Fn := (mem_sdiff.mp hj_in_diff).1
  use j, hj_in_Fn
  intro x hx
  have hx_Fn : x ∈ Fn := (mem_sdiff.mp hx).1
  have hx_ne_j : x ∉ ({j} : Finset α) := (mem_sdiff.mp hx).2
  have hx_not_diff : x ∉ Fn \ (Fn ∩ F_star) := by
    rw [hj]
    exact hx_ne_j
  have hx_inter : x ∈ Fn ∩ F_star := by
    by_contra hc
    apply hx_not_diff
    exact mem_sdiff.mpr ⟨hx_Fn, hc⟩
  exact (mem_inter.mp hx_inter).2

/-- Proposition 1: Cut validity with a Global Lower Bound (GLB).
    If L is a valid global lower bound (L ≤ ν(F*) for all F*),
    and adding jobs never decreases setup costs (subproblem monotonicity under superset),
    then cut (27) is valid for ANY master problem solution F*. -/
theorem cut_validity_with_GLB
    (Fn : Finset α)
    (ν : Finset α → ℝ)
    (θ_hat L : ℝ)
    (h_pos : 0 ≤ θ_hat - L)
    (h_superset : ∀ F_star : Finset α, Fn ⊆ F_star → θ_hat ≤ ν F_star)
    (h_glb : ∀ F_star : Finset α, L ≤ ν F_star)
    (F_star : Finset α) :
    cutRHS Fn.card θ_hat L (Fn ∩ F_star).card ≤ ν F_star := by
  have hκ_le : (Fn ∩ F_star).card ≤ Fn.card := card_le_card inter_subset_left
  by_cases h_full : (Fn ∩ F_star).card = Fn.card
  · -- Case 1: κ = |Fn|, which means Fn ⊆ F_star
    have h_sub : Fn ⊆ F_star := by
      have h_eq : Fn ∩ F_star = Fn := eq_of_subset_of_card_le inter_subset_left (by omega)
      intro x hx
      have : x ∈ Fn ∩ F_star := by rw [h_eq]; exact hx
      exact (mem_inter.mp this).2
    rw [h_full, cutRHS_at_card]
    exact h_superset F_star h_sub
  · -- Case 2: κ ≤ |Fn| - 1
    have hκ_lt : (Fn ∩ F_star).card + 1 ≤ Fn.card := by omega
    rw [cutRHS_eq]
    have h_coeff : ((Fn ∩ F_star).card : ℝ) - (Fn.card : ℝ) + 1 ≤ 0 := by
      have : ((Fn ∩ F_star).card : ℝ) + 1 ≤ (Fn.card : ℝ) := by exact_mod_cast hκ_lt
      linarith
    have h_neg : (θ_hat - L) * (((Fn ∩ F_star).card : ℝ) - (Fn.card : ℝ) + 1) ≤ 0 := by
      nlinarith
    have h_rhs_le_L : (θ_hat - L) * (((Fn ∩ F_star).card : ℝ) - (Fn.card : ℝ) + 1) + L ≤ L := by
      linarith
    exact le_trans h_rhs_le_L (h_glb F_star)

/-!
### Proof of Proposition 2: Local Lower Bound (LLB)
-/

/-- Step 4 Algebra:
    For any m ≥ 1 and ψ(Fn) ≤ θ_hat:
    θ_hat - m * (θ_hat - LLB) ≤ ψ(Fn) - m * (ψ(Fn) - LLB).
    The difference is (m - 1) * (θ_hat - ψ(Fn)) ≥ 0. -/
lemma step4_bound (ψ_n θ_hat LLB : ℝ) (m : ℕ) (hm : 1 ≤ m) (hθ : ψ_n ≤ θ_hat) :
    θ_hat - (m : ℝ) * (θ_hat - LLB) ≤ ψ_n - (m : ℝ) * (ψ_n - LLB) := by
  have h_diff : (ψ_n - (m : ℝ) * (ψ_n - LLB)) - (θ_hat - (m : ℝ) * (θ_hat - LLB)) =
                ((m : ℝ) - 1) * (θ_hat - ψ_n) := by ring
  have hm_ge : 0 ≤ (m : ℝ) - 1 := by
    have : 1 ≤ (m : ℝ) := by exact_mod_cast hm
    linarith
  have hθ_ge : 0 ≤ θ_hat - ψ_n := by linarith
  have h_prod : 0 ≤ ((m : ℝ) - 1) * (θ_hat - ψ_n) := mul_nonneg hm_ge hθ_ge
  linarith

/-- Step 3 Induction:
    Bounding the decrease of ψ when a subset R of m jobs is removed from Fn.
    ψ(Fn \ R) ≥ ψ(Fn) - m * (ψ(Fn) - LLB). -/
lemma step3_induction
    (Fn : Finset α)
    (ψ : Finset α → ℝ)
    (LLB : ℝ)
    (h_LLB : ∀ j ∈ Fn, LLB ≤ ψ (Fn \ {j}))
    (h_marginal : ∀ F' ⊆ Fn, ∀ j ∈ F', ψ F' - ψ (F' \ {j}) ≤ ψ Fn - ψ (Fn \ {j}))
    (R : Finset α) (hR : R ⊆ Fn) :
    ψ Fn - (R.card : ℝ) * (ψ Fn - LLB) ≤ ψ (Fn \ R) := by
  induction R using Finset.induction_on with
  | empty =>
    simp
  | @insert j R' hj ih =>
    have hj_in_Fn : j ∈ Fn := hR (mem_insert_self j R')
    have hR'_sub : R' ⊆ Fn := (subset_insert j R').trans hR
    have ih' := ih hR'_sub
    have h_diff_eq : Fn \ (insert j R') = (Fn \ R') \ {j} := by
      ext x
      simp only [mem_sdiff, mem_insert, mem_singleton, not_or]
      tauto
    rw [h_diff_eq]
    have hF'_sub : Fn \ R' ⊆ Fn := sdiff_subset
    have hj_in_F' : j ∈ Fn \ R' := mem_sdiff.mpr ⟨hj_in_Fn, hj⟩
    have h_step := h_marginal (Fn \ R') hF'_sub j hj_in_F'
    have h_llb_step : ψ Fn - ψ (Fn \ {j}) ≤ ψ Fn - LLB := by
      have := h_LLB j hj_in_Fn
      linarith
    have h_bound_step : ψ (Fn \ R') - ψ ((Fn \ R') \ {j}) ≤ ψ Fn - LLB := by
      linarith
    have h_card : ((insert j R').card : ℝ) = (R'.card : ℝ) + 1 := by
      rw [card_insert_of_notMem hj]
      push_cast
      ring
    rw [h_card]
    linarith

/-- Master Theorem: Formal Verification of Proposition 2.
    Under:
    - Step 1: Valid subproblem lower bounding surrogate ψ (ψ G ≤ ν G),
    - Step 2: Monotonicity under set inclusion / metric shortcutting (A ⊆ B → ν A ≤ ν B),
    - Optimality: ψ Fn ≤ θ_hat and θ_hat ≤ ν Fn,
    - Definition of LLB: LLB ≤ ψ(Fn \ {j}) for all j ∈ Fn,
    - Supermodularity / marginal decrease property of ψ (Step 3),
    the lifted cut (27) with lower bound LLB is valid for ALL master problem solutions F*. -/
theorem proposition2_cut_validity
    (Fn : Finset α)
    (ν : Finset α → ℝ)
    (ψ : Finset α → ℝ)
    (θ_hat LLB : ℝ)
    (h_mono : ∀ A B : Finset α, A ⊆ B → ν A ≤ ν B)
    (h_bound : ∀ G : Finset α, ψ G ≤ ν G)
    (h_opt : ψ Fn ≤ θ_hat)
    (h_subprob_opt : θ_hat ≤ ν Fn)
    (h_LLB : ∀ j ∈ Fn, LLB ≤ ψ (Fn \ {j}))
    (h_marginal : ∀ F' ⊆ Fn, ∀ j ∈ F', ψ F' - ψ (F' \ {j}) ≤ ψ Fn - ψ (Fn \ {j}))
    (F_star : Finset α) :
    cutRHS Fn.card θ_hat LLB (Fn ∩ F_star).card ≤ ν F_star := by
  let G := Fn ∩ F_star
  let R := Fn \ G
  have hG_sub_Fn : G ⊆ Fn := inter_subset_left
  have hG_sub_F_star : G ⊆ F_star := inter_subset_right
  have hR_sub_Fn : R ⊆ Fn := sdiff_subset
  have hFn_sdiff_R : Fn \ R = G := Finset.sdiff_sdiff_eq_self hG_sub_Fn
  have hR_card : R.card = Fn.card - G.card := card_sdiff_of_subset hG_sub_Fn
  have hG_card_le : G.card ≤ Fn.card := card_le_card hG_sub_Fn
  have h_rhs_eq : cutRHS Fn.card θ_hat LLB G.card =
      θ_hat - (R.card : ℝ) * (θ_hat - LLB) := by
    dsimp [cutRHS]
    rw [hR_card]
    have h_cast : ((Fn.card - G.card : ℕ) : ℝ) = (Fn.card : ℝ) - (G.card : ℝ) :=
      Nat.cast_sub hG_card_le
    rw [h_cast]
    ring
  rw [h_rhs_eq]
  by_cases hm0 : R.card = 0
  · -- Case m = 0: F_star contains Fn as a subset
    have hG_eq_Fn : G = Fn := eq_of_subset_of_card_le hG_sub_Fn (by omega)
    have hFn_sub : Fn ⊆ F_star := by
      intro x hx
      have : x ∈ G := by rw [hG_eq_Fn]; exact hx
      exact (mem_inter.mp this).2
    rw [hm0]
    push_cast
    ring_nf
    have : ν Fn ≤ ν F_star := h_mono Fn F_star hFn_sub
    linarith
  · -- Case m ≥ 1: Chaining Steps 2, 1, 3, and 4
    have hm1 : 1 ≤ R.card := by omega
    have h_step4 := step4_bound (ψ Fn) θ_hat LLB R.card hm1 h_opt
    have h_step3 := step3_induction Fn ψ LLB h_LLB h_marginal R hR_sub_Fn
    rw [hFn_sdiff_R] at h_step3
    have h_step1 : ψ G ≤ ν G := h_bound G
    have h_step2 : ν G ≤ ν F_star := h_mono G F_star hG_sub_F_star
    linarith
