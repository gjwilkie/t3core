#!/bin/julia
using species: init_species
using input: Nrad, Nv, read_input, deltat,Nt, dilution_model, diffmodel, maxwellian_edge
using matrix: init_collop, build_matrix, solve_steadystate, f0, gindex, advance_timestep
using diffcoeff: init_diffCoeff, init_diffCoeff_zero
using sourcemod: init_alpha_source
using boundary: calculate_boundary, F0edge
using geometry: init_geometry
using grids: init_rgrid, v, init_vgrid
using postproc: plot_steadystate, init_postproc, save_f0, save_density

function main()

if length(ARGS) == 1
  runname = ascii(ARGS[1])
elseif length(ARGS) == 0
  runname = "test"
else
  error("Cannot parse more than one argument")
end
try run(`mkdir $runname`)
catch
  println("Directory "*runname*"already exists")
end
  
println(runname)

read_input()

# Create radial grid
init_rgrid()

# From GS2 netcdf files, read in species densities, temperatures, etc.
init_species()

# Read geometric properties from GS2 files
init_vgrid()

# Read geometric properties from GS2 files
init_geometry()

# Initialize differentiation weights
init_collop() 

# Initialize source
init_alpha_source() 

# Calculate the diffusion coefficients from GS2 files
init_diffCoeff()

# Construct global Matrix
build_matrix()

#if !semianalytic_on
  # Calculate F0 at the boundary
#end
calculate_boundary()

init_postproc(runname)

if deltat>0.0
  for it in 1:Nt
    println("Step ",it)
    advance_timestep(it)
    save_f0(it)
  end 
else
  # Solve matrix equation
  solve_steadystate()
  save_f0(1)
end

if dilution_model > 0
   save_density()
end
if deltat < 0.0
  plot_steadystate()
end

end

main()

