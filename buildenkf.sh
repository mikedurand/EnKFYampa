#!/bin/sh -x

# Useful Switches:
# -C    : enables all checks on run-time conditions
# -CB   : enables checks on array bounds only
# -g    : produces symbolic debug information
# -fpe0 : Floating point underflow results in zero; all others abort execution
# -traceback: Traces where the error occurs: Yippee!
# -fpstkchk : Facilitates locating the source of floating point exceptions

mpif90 enkf.f90 run_params.f90 rsvd.f90 mvnrnd.f90 rand_normal.f90 ssib3.f compileystats.f90 interfacez.f90 ssib_layer.f90 rad_xfer_model.f90 atm_model.f90 can_model.f90 ss_model.f90 kupdate.f90 invert_matrix.f90 serial_date.f90 std2j.f90 leapyear.f90 disaggregate.f90 get_min_loc.f90 read_veg.f90 vegin.f get_sla_data.f90 get_veg_data.f cholesky.f90 interpolate.f90 -o ./filter
# rm *.o
