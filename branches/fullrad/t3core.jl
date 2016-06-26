#!/bin/julia
using species: init_species
using input: Nrad, Nv, read_input, deltat,Nt, dilution_model, diffmodel, maxwellian_edge, Nout
using matrix: init_collop, build_matrix, solve_steadystate, f0, gindex, advance_timestep
using diffcoeff: init_diffCoeff, init_diffCoeff_zero
using sourcemod: init_alpha_source
using boundary: calculate_boundary, F0edge
using geometry: init_geometry
using grids: init_rgrid, v, init_vgrid
using postproc: plot_steadystate, init_postproc, save_density, transient_diagnostics

function main()

verbose = false

if length(ARGS) == 1
  runname = ascii(ARGS[1])
elseif length(ARGS) == 0
  runname = "test"
else
  error("Cannot parse more than one argument")
end
try run(`mkdir $runname`)
catch
end
  
read_input()

# Create radial grid
verbose && println("Creating radial grid...")
init_rgrid()

# From GS2 netcdf files, read in species densities, temperatures, etc.
verbose && println("Setting up species data...")
init_species()

# Read geometric properties from GS2 files
verbose && println("Creating velocity grid...")
init_vgrid()

# Read geometric properties from GS2 files
verbose && println("Setting up geometry...")
init_geometry()

# Initialize differentiation weights
verbose && println("Calculating collision operator...")
init_collop() 

# Initialize source
verbose && println("Calculating source...")
init_alpha_source() 

# Calculate the diffusion coefficients from GS2 files
verbose && println("Calculating diffusion coefficients...")
init_diffCoeff()

# Construct global Matrix
verbose && println("Assembling matrix...")
build_matrix()

#if !semianalytic_on
  # Calculate F0 at the boundary
#end
verbose && println("Calculating edge distribution...")
calculate_boundary()

verbose && println("Outputting diagnostics...")
init_postproc(runname)

if deltat>0.0
  for it in 1:Nt
    println("Step ",it)
    advance_timestep(it)
    if in(it,unique(round(linspace(1,Nt,Nout))))
      transient_diagnostics(f0,it)
    end
  end 
else
  # Solve matrix equation
  verbose && println("Inverting matrix...")
  solve_steadystate()
end

if dilution_model > 0
   save_density(f0)
end
verbose && println("Processing data...")
if deltat < 0.0
  plot_steadystate(f0)
end

end

main()

