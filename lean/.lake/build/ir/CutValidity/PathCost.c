// Lean compiler output
// Module: CutValidity.PathCost
// Imports: public import Init public meta import Init public import Mathlib.Data.List.Basic public import Mathlib.Tactic
#include <lean/lean.h>
#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-label"
#elif defined(__GNUC__) && !defined(__CLANG__)
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-label"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#ifdef __cplusplus
extern "C" {
#endif
extern lean_object* lp_mathlib___private_Mathlib_Basic_Real_Basic_0__Real_zero;
lean_object* lp_mathlib___private_Mathlib_Basic_Real_Basic_0__Real_add___lam__0(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_CutValidity_pathCost___redArg(lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_CutValidity_pathCost(lean_object*, lean_object*, lean_object*);
LEAN_EXPORT lean_object* lp_CutValidity_pathCost___redArg(lean_object* v_d_1_, lean_object* v_x_2_){
_start:
{
if (lean_obj_tag(v_x_2_) == 0)
{
lean_object* v___x_3_; 
lean_dec(v_d_1_);
v___x_3_ = lp_mathlib___private_Mathlib_Basic_Real_Basic_0__Real_zero;
return v___x_3_;
}
else
{
lean_object* v_tail_4_; 
v_tail_4_ = lean_ctor_get(v_x_2_, 1);
lean_inc(v_tail_4_);
if (lean_obj_tag(v_tail_4_) == 0)
{
lean_object* v___x_5_; 
lean_dec_ref_known(v_x_2_, 2);
lean_dec(v_d_1_);
v___x_5_ = lp_mathlib___private_Mathlib_Basic_Real_Basic_0__Real_zero;
return v___x_5_;
}
else
{
lean_object* v_head_6_; lean_object* v_head_7_; lean_object* v___x_8_; lean_object* v___x_9_; lean_object* v___f_10_; 
v_head_6_ = lean_ctor_get(v_x_2_, 0);
lean_inc(v_head_6_);
lean_dec_ref_known(v_x_2_, 2);
v_head_7_ = lean_ctor_get(v_tail_4_, 0);
lean_inc(v_d_1_);
lean_inc(v_head_7_);
v___x_8_ = lean_apply_2(v_d_1_, v_head_6_, v_head_7_);
v___x_9_ = lp_CutValidity_pathCost___redArg(v_d_1_, v_tail_4_);
v___f_10_ = lean_alloc_closure((void*)(lp_mathlib___private_Mathlib_Basic_Real_Basic_0__Real_add___lam__0), 3, 2);
lean_closure_set(v___f_10_, 0, v___x_8_);
lean_closure_set(v___f_10_, 1, v___x_9_);
return v___f_10_;
}
}
}
}
LEAN_EXPORT lean_object* lp_CutValidity_pathCost(lean_object* v_V_11_, lean_object* v_d_12_, lean_object* v_x_13_){
_start:
{
lean_object* v___x_14_; 
v___x_14_ = lp_CutValidity_pathCost___redArg(v_d_12_, v_x_13_);
return v___x_14_;
}
}
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_Init(uint8_t builtin);
lean_object* initialize_mathlib_Mathlib_Data_List_Basic(uint8_t builtin);
lean_object* initialize_mathlib_Mathlib_Tactic(uint8_t builtin);
void lean_initialize();
static bool _G_initialized = false;
LEAN_EXPORT lean_object* initialize_CutValidity_CutValidity_PathCost(uint8_t builtin) {
lean_object * res;
if (_G_initialized) return lean_io_result_mk_ok(lean_box(0));
_G_initialized = true;
lean_initialize();
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_Init(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_mathlib_Mathlib_Data_List_Basic(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
res = initialize_mathlib_Mathlib_Tactic(builtin);
if (lean_io_result_is_error(res)) return res;
lean_dec_ref(res);
return lean_io_result_mk_ok(lean_box(0));
}
#ifdef __cplusplus
}
#endif
