SETUP_FLAGS := -d=22
SCHEME_DIR := ../../schemes
TESTS := impl_2step_tvdlf_mm.log impl_3step_hll_cada.log impl_4step_hll_mc.log	\
impl_4step_hllc_ko.log impl_rk4_tvdlf_cada.log

impl_2step_tvdlf_mm.log: test.par $(SCHEME_DIR)/2step_tvdlf_mm.par
impl_3step_hll_cada.log: test.par $(SCHEME_DIR)/3step_hll_cada.par
impl_4step_hll_mc.log: test.par $(SCHEME_DIR)/4step_hll_mc.par
impl_4step_hllc_ko.log: test.par $(SCHEME_DIR)/4step_hllc_ko.par
impl_rk4_tvdlf_cada.log: test.par $(SCHEME_DIR)/rk4_tvdlf_cada.par

include ../../test_rules.make
