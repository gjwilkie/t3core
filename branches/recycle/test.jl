using input: input_for_tests, Nrad, Nv, vmax, input_for_analytic_test
using grids: init_rgrid, init_vgrid, v, ddv
using species: init_species_for_tests
using geometry: init_geometry
using matrix: init_collop, build_matrix, collop, zero_collop, gindex, global_matrix, Vprime_global
using sourcemod: init_alpha_source, source_local
using boundary: calculate_boundary, F0edge
using diffcoeff: init_diffCoeff_zero
using Base.Test
using Winston
using constants: valpha

# Do analytic test:
input_for_tests()

print("Testing input... ")
@test_approx_eq 15474962.7102508 vmax
print("Success!\n")

init_rgrid()

init_species_for_tests()

print("Testing v grid... ")
init_vgrid( )
@test_approx_eq 7494441.19240807 v[99]
@test_approx_eq 1.27869276765246e-5 ddv[45,46]
print("Success!\n")

init_geometry()

print("Testing collision operator... ") 
init_collop()
@test_approx_eq 128.088688125756 collop[50,30,1]
@test_approx_eq -.498913649121981 collop[135,73,1]
print("Success!\n")

print("Testing source... ") 
init_alpha_source()
@test_approx_eq 4.74764299169435e-4 source_local[end,150]
print("Success!\n")

init_diffCoeff_zero()

print("Testing matrix... ")
build_matrix()
idx = gindex(135,1)
jdx = gindex(73,1)
@test_approx_eq .498913649121981*Vprime_global[73,73] global_matrix[idx,jdx]
print("Success!\n")

print("Testing F0edge... ")
calculate_boundary()
@test_approx_eq_eps 1.7879918175437e-4 F0edge[70] 1.0e-10
@test_approx_eq_eps 1.60720325835966e-5 F0edge[134] 1.0e-10
print("Success!\n")


#semilogy(v/valpha,abs(F0edge))

# Be sure to test:
# - f0edge
# - Collision operator
# - Analytic coupled r-v solution (w/o collisions)
