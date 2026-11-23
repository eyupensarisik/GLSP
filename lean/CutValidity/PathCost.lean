import Mathlib.Data.List.Basic
import Mathlib.Tactic

set_option linter.style.header false
set_option linter.style.longLine false
set_option linter.unusedSectionVars false

/-!
# Metric Path Costs and Local Lower Bound Verification

This file provides foundational results for sequencing problems with
sequence-dependent setup costs satisfying the triangle inequality:
1. Metric path cost definitions.
2. Shortcutting lemmas: vertex removal never increases path cost.
3. Formal proof of Equation (29): within-period lower bound based on minimum outgoing setup costs.
-/

variable {V : Type*} (d : V → V → ℝ)

/-- The total setup cost of processing a sequence `l` of jobs on a machine. -/
def pathCost : List V → ℝ
  | [] => 0
  | [_] => 0
  | u :: v :: rest => d u v + pathCost (v :: rest)

lemma pathCost_nil : pathCost d [] = 0 := rfl
lemma pathCost_singleton (u : V) : pathCost d [u] = 0 := rfl
lemma pathCost_cons_cons (u v : V) (rest : List V) :
    pathCost d (u :: v :: rest) = d u v + pathCost d (v :: rest) := rfl

/-- Path cost is non-negative whenever the distance matrix is non-negative. -/
lemma pathCost_nonneg (h_nonneg : ∀ u v, 0 ≤ d u v) (l : List V) :
    0 ≤ pathCost d l := by
  induction l with
  | nil => simp [pathCost]
  | cons u rest ih =>
    cases rest with
    | nil => simp [pathCost]
    | cons v r =>
      change 0 ≤ d u v + pathCost d (v :: r)
      have := h_nonneg u v
      have : 0 ≤ pathCost d (v :: r) := ih
      linarith

/-- Dropping the first job from a sequence never increases path cost. -/
lemma pathCost_tail (h_nonneg : ∀ u v, 0 ≤ d u v) (v : V) (tail : List V) :
    pathCost d tail ≤ pathCost d (v :: tail) := by
  cases tail with
  | nil => simp [pathCost]
  | cons w rest =>
    simp [pathCost]
    have := h_nonneg v w
    linarith

/-- Dropping the second job from a sequence never increases path cost under triangle inequality. -/
lemma pathCost_skip_second (h_nonneg : ∀ u v, 0 ≤ d u v)
    (h_tri : ∀ u v w, d u w ≤ d u v + d v w) (u v : V) (tail : List V) :
    pathCost d (u :: tail) ≤ pathCost d (u :: v :: tail) := by
  cases tail with
  | nil =>
    simp [pathCost]
    have := h_nonneg u v
    linarith
  | cons w rest =>
    simp [pathCost]
    have := h_tri u v w
    linarith

/-- Metric Shortcutting Lemma:
    Removing any job `x` at any position from a sequence `l₁ ++ [x] ++ l₂`
    never increases the total setup cost. -/
theorem pathCost_shortcut (h_nonneg : ∀ u v, 0 ≤ d u v)
    (h_tri : ∀ u v w, d u w ≤ d u v + d v w)
    (l₁ l₂ : List V) (x : V) :
    pathCost d (l₁ ++ l₂) ≤ pathCost d (l₁ ++ [x] ++ l₂) := by
  induction l₁ with
  | nil =>
    dsimp
    exact pathCost_tail d h_nonneg x l₂
  | cons a tail ih =>
    cases tail with
    | nil =>
      dsimp
      exact pathCost_skip_second d h_nonneg h_tri a x l₂
    | cons b rest =>
      change d a b + pathCost d (b :: rest ++ l₂) ≤ d a b + pathCost d (b :: rest ++ [x] ++ l₂)
      have h_ih := ih
      linarith

variable (min_sc : V → ℝ)

/-- In any sequence of distinct jobs, internal transitions pay at least the sum of
    outgoing minimum setup costs for all jobs except the last one. -/
lemma pathCost_ge_sum_min_sc_dropLast
    (h_min : ∀ u v, u ≠ v → min_sc u ≤ d u v) :
    ∀ l : List V, l.Nodup →
    ((l.dropLast).map min_sc).sum ≤ pathCost d l
  | [] => by simp [pathCost]
  | [u] => by simp [pathCost]
  | u :: v :: rest => by
    intro h_nodup
    have ⟨h_not_mem, h_rest_nodup⟩ := List.nodup_cons.mp h_nodup
    have h_ne : u ≠ v := by
      intro heq
      apply h_not_mem
      simp [heq]
    have ih := pathCost_ge_sum_min_sc_dropLast h_min (v :: rest) h_rest_nodup
    have h_drop : (u :: v :: rest).dropLast = u :: (v :: rest).dropLast := rfl
    rw [h_drop]
    simp only [List.map_cons, List.sum_cons, pathCost]
    have h_step : min_sc u ≤ d u v := h_min u v h_ne
    linarith

/-- Formal Verification of Equation (29):
    For any sequence `l` of distinct products in a period, the setup cost is bounded below
    by the sum of minimum setup costs minus the maximum minimum setup cost in `l`. -/
theorem pathCost_ge_sum_sub_max
    (h_min : ∀ u v, u ≠ v → min_sc u ≤ d u v)
    (l : List V) (h_nodup : l.Nodup) (h_nonempty : l ≠ [])
    (M : ℝ) (hM : ∀ x ∈ l, min_sc x ≤ M) :
    (l.map min_sc).sum - M ≤ pathCost d l := by
  have h_bound := pathCost_ge_sum_min_sc_dropLast d min_sc h_min l h_nodup
  have h_sum_eq : (l.map min_sc).sum = ((l.dropLast).map min_sc).sum + min_sc (l.getLast h_nonempty) := by
    conv_lhs => rw [← List.dropLast_append_getLast h_nonempty]
    simp [List.map_append]
  have h_last_le_M : min_sc (l.getLast h_nonempty) ≤ M := by
    apply hM
    exact List.getLast_mem h_nonempty
  linarith
