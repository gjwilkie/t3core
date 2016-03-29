#!/bin/julia
using species: init_species
using input: Nrad, Nv, read_input, deltat,Nt, dilution_model, semianalytic_on, a
using matrix: init_collop, build_matrix, solve_steadystate, f0, gindex, advance_timestep, build_matrix_maxw
using diffcoeff: init_diffCoeff, init_diffCoeff_zero, Dnn, DnT, DTn, DTT, init_diffCoeff_maxw, pflux0, hflux0
using sourcemod: init_alpha_source, init_maxw_source, source_local
using boundary: calculate_boundary, maxwellian_boundary, F0edge, calculate_boundary_maxw
using geometry: init_geometry
using grids: init_rgrid, v, init_vgrid, rgrid
using postproc: plot_f0, plot_descriminant, init_postproc, save_f0, plot_alphaprofile, plot_source, save_density, plot_fluxarrows
using PyPlot

plt[:style][:use]("myfig")

read_input()

# Create radial grid
println("Creating radial grid...")
init_rgrid()

# From GS2 netcdf files, read in species densities, temperatures, etc.
println("Setting up species data...")
init_species()

# Read geometric properties from GS2 files
println("Creating velocity grid...")
init_vgrid()

# Read geometric properties from GS2 files
println("Setting up geometry...")
init_geometry()

# Initialize source
println("Calculating source...")
init_maxw_source() 

# Calculate the diffusion coefficients from GS2 files
println("Calculating diffusion coefficients...")
if Nrad>1 
  init_diffCoeff_maxw()
end

plot(rgrid/a,Dnn)
xlabel(L"$r / a$")
savefig("Dnn.png")
cla()
clf()
plot(rgrid/a,DnT)
xlabel(L"$r / a$")
savefig("DnT.png")
cla()
clf()
plot(rgrid/a,DTn)
xlabel(L"$r / a$")
savefig("DTn.png")
cla()
clf()
plot(rgrid/a,DTT)
xlabel(L"$r / a$")
savefig("DTT.png")
cla()
clf()
plot(rgrid/a,pflux0)
xlabel(L"$r / a$")
savefig("pflux0.png")
cla()
clf()
plot(rgrid/a,hflux0)
xlabel(L"$r / a$")
savefig("hflux0.png")
cla()
clf()

plot(rgrid/a,source_local)
xlabel(L"$r / a$")
savefig("source.png")
cla()
clf()

# Construct global Matrix
println("Assembling matrix...")
build_matrix_maxw()

println("Setting up boundary condition...")
calculate_boundary_maxw()

solve_steadystate()

density_maxw = f0[1:Nrad]
temp_maxw = f0[Nrad+1:2*Nrad]

plot(rgrid/a,density_maxw)
xlabel(L"$r / a$")
ylabel(L"$n$")
savefig("density.png")
cla()
clf()

plot(rgrid/a,temp_maxw)
xlabel(L"$r / a$")
ylabel(L"$T$")
savefig("temperature.png")
cla()
clf()
