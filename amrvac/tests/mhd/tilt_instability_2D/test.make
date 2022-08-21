SETUP_FLAGS := -d=2
SCHEME_DIR := ../../schemes
SCHEMES := 2step_tvdlf_mm 2step_tvdmu_al 3step_hll_cada 4step_hll_mc \
4step_hllc_ko rk4_tvdlf_cada ssprk54_fd_mp5 ssprk54_hll_mp5 3step_hlld_cada \
ssprk54_hlld_mp5

TESTS := $(SCHEMES:%=ti_2d_%.log)

include ../../test_rules.make

# Generate dependency rules for the tests
$(foreach s, $(SCHEMES),\
	$(eval $(s:%=ti_2d_%.log): ti_2d.par $(SCHEME_DIR)/$(s).par))
