module species
using NetCDF
using input: Nrad_gs2, mref, qref, tracespecs, Nrad, rgrid_gs2, Te_in, Ti_in, ne_in, rgrid_in, DTmix, m_trace, Z_trace, spline_k
using Dierckx
using grids: rgrid,v,d3v,ddv
using constants: me, el, mp

export init_species,nspec,charge,mass,bulkspecs, dens_gs2, temp_gs2, filenames, vcrit, ne

filenames=String[]
vcrit=Float64[]
mass=Float64[]
charge=Float64[]
dens_gs2 =Float64[]
temp_gs2 =Float64[]
ne =Float64[]
Te =Float64[]
Ti =Float64[]
nref =Float64[]
Tref =Float64[]
bulkspecs=Int64[]

function init_species()
  # Create strings for filenames
  # files are located in nc directory and are named:
  # nc/r1.nc, r2.nc, r3.nc, etc.
  global filenames, nspec, bulkspecs, charge, mass, vcrit, temp_gs2, dens_gs2, ne, Te, Ti, nref, Tref

  for i in 1:Nrad_gs2
    push!(filenames,"nc/r$i.nc") 
  end
 
  # Find and check number of species
  nspec = ncread(filenames[1],"nspecies")[1]
  for i in 1:Nrad_gs2
    if (ncread(filenames[i],"nspecies")[1] != nspec)
       error("All files must have the same number of bulk species")
    end
  end
  
  # Account bulk species separately
  for ispec in 1:nspec
    if !(ispec in tracespecs)
      push!(bulkspecs,ispec)
    end
  end  

  # Read species parameters
  mass=ncread(filenames[1],"mass")[:]*mref
  charge=ncread(filenames[1],"charge")[:]*qref

  # Radial-dependent equilibrium species quantities
  dens_gs2 = zeros(Nrad_gs2,nspec)
  temp_gs2 = zeros(Nrad_gs2,nspec)
  nref = zeros(Nrad_gs2)
  Tref = zeros(Nrad_gs2)
  vcrit = zeros(Nrad)

  for is in 1:nspec
    # Find the density and temperature profiles
    nref_func = Spline1D(rgrid_in, ne_in,k=spline_k)
    Tref_func = Spline1D(rgrid_in, Ti_in,k=spline_k)
    for ir in 1:Nrad_gs2
      nref[ir] = evaluate(nref_func,rgrid_gs2[ir]) 
      Tref[ir] = evaluate(Tref_func,rgrid_gs2[ir]) 
      dens_gs2[ir,is]=ncread(filenames[ir],"dens")[is]*nref[ir]
      temp_gs2[ir,is]=ncread(filenames[ir],"temp")[is]*Tref[ir]
    end

  end

  ne = zeros(Float64,Nrad)
  Te = zeros(Float64,Nrad)
  Ti = zeros(Float64,Nrad)

  ne_func = Spline1D(rgrid_in,ne_in,k=spline_k)
  Te_func = Spline1D(rgrid_in,Te_in,k=spline_k)
  Ti_func = Spline1D(rgrid_in,Ti_in,k=spline_k)
  for ir in 1:Nrad
    ne[ir] = evaluate(ne_func,rgrid[ir])
    Te[ir] = evaluate(Te_func,rgrid[ir])
    Ti[ir] = evaluate(Ti_func,rgrid[ir])
  end

  for ir in 1:Nrad
    ZI = DTmix*ne[ir]*m_trace/(ne[ir]*2.0*mp)
    ZI += (1.0-DTmix)*ne[ir]*m_trace/(ne[ir]*3.0*mp)

    vcrit[ir] = (3.0*sqrt(pi)*me*ZI/(4.0*m_trace))^(1/3)*sqrt(2.0*Te[ir]/me)
  end
 
end

function init_species_for_analytic_test()

end

end
