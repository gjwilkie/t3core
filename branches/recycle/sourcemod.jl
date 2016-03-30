module sourcemod
using input: Nrad, Nv, tracespecs, a, diffmodel, DTmix, rgrid_in, Ti_in, ne_in, m_trace, D0, brv,bvr, rmaj, zerosource, dilution_model, ashmode, Tashfac
using grids: v, d3v, rgrid, v
using constants: mp, ev2joule, Ealpha, el, valpha
using species: ne, Ti, mass
using geometry: Vprime
using Dierckx

export init_alpha_source, source, add2source, source_local, zero_source_element, reaction_rate, add2source_element

Tavg = Float64[] 
source = Float64[]
source_local = Float64[]
source_in = Float64[]
reaction_rate = Float64[]

# Source for alpha particles
function init_alpha_source()
  global source, source_local, Tavg, reaction_rate, source_in

  source= zeros(Nrad*Nv)
  source_local = zeros(Nrad,Nv)
  source_in = zeros(length(rgrid_in),Nv)
  reaction_rate = zeros(Nrad)
  reaction_rate_in = zeros(length(rgrid_in))

  nhat = Array(Float64,Nrad)
  if ( (dilution_model == 2 ) || (dilution_model == 3)) & isfile("density.dat")
    data = readdlm("density.dat",' ')
    rgrid_d = data[:,1]
    dens_d = data[:,2]
    dens_func = Spline1D(rgrid_d,dens_d)
    nhat = evaluate(dens_func,rgrid)
  else
    nhat[:] = 0.0
  end

  for ir in 1:length(rgrid_in)
    # Express Ti in kev
    Ti_kev = Ti_in[ir]/(1000.0*el)
   
    # If diluting, take away some density of deuterium and tritium (half of alpha density each)
    ni2 = (DTmix*(1.0-nhat[ir])ne_in[ir])*((1.0-DTmix)*(1.0-nhat[ir])*ne_in[ir])
  
    # Find fusion reaction rate
    # From the NRL plasma formulary

    reaction_rate_in[ir] = ni2 * 3.68e-12 * Ti_kev^(-2.0/3.0) * exp(-19.94 * Ti_kev^(-1.0/3.0)) * 1.0e-6
     
    integrated_source = 0.0
    for iv in 1:Nv
      b = 5.0 / 16.0
      E = 0.5 * m_trace * v[iv]^2 
      if ashmode
        source_in[ir,iv] = exp(-m_trace*v[iv]^2/(2.0*Ti_in[ir]*Tashfac))
      else
        source_in[ir,iv] = exp(-b*(E - Ealpha)^2/(Ti_in[ir]*Ealpha))
      end
      integrated_source += source_in[ir,iv]*d3v[iv]
    end

    # Normalize the source appropriately
    for iv in 1:Nv
      source_in[ir,iv] = source_in[ir,iv]*reaction_rate_in[ir]/integrated_source
    end
    
  end


  for ir in 1:Nrad
    # Express Ti in kev
    Ti_kev = Ti[ir]/(1000.0*el)
   
    ni2 = (DTmix*(1.0-nhat[ir])*ne[ir])*((1.0-DTmix)*(1.0-nhat[ir])*ne[ir])
  
    # Find fusion reaction rate
    # From the NRL plasma formulary

    reaction_rate[ir] = ni2 * 3.68e-12 * Ti_kev^(-2.0/3.0) * exp(-19.94 * Ti_kev^(-1.0/3.0)) * 1.0e-6
     
    integrated_source = 0.0
    for iv in 1:Nv
      idx = gindex(ir,iv)
      b = 5.0 / 16.0
      E = 0.5 * m_trace * v[iv]^2 
      if ashmode
        source[idx] = exp(-m_trace*v[iv]^2/(2.0*Ti[ir]*Tashfac))
      else
        source[idx] = exp(-b*(E - Ealpha)^2/(Ti[ir]*Ealpha))
      end
      integrated_source += source[idx]*d3v[iv]
    end

    # Normalize the source appropriately
    for iv in 1:Nv
      idx = gindex(ir,iv)
      source[idx] = source[idx]*reaction_rate[ir]/integrated_source
      source_local[ir,iv] = source[idx]  #< This array is used to calculate the SD distribution at a given radius, and does not include the V'(psi)*v^2 factor
      source[idx] *= Vprime[ir]*v[iv]^2
    end

  end
  
  for iv in 1:Nv
     idx = gindex(Nrad,iv)
     source[idx] = 0.0
     source_local[Nrad,iv] = 0.0
  end

  if zerosource
     source[:] = 0.0
     source_local[:,:] = 0.0
  end

end

function init_maxw_source()
  global source, source_local, Tavg, reaction_rate, source_in

  source= zeros(Nrad*2)
  source_local = zeros(Nrad)
  source_in = zeros(length(rgrid_in))
  reaction_rate = zeros(Nrad)
  reaction_rate_in = zeros(length(rgrid_in))

  # First define the source at all (input) radii. Used to calculate boundary condition.
  for ir in 1:length(rgrid_in)
    # Express Ti in kev
    Ti_kev = Ti_in[ir]/(1000.0*ev2joule)
   
    ni2 = (DTmix*ne_in[ir])*((1.0-DTmix)*ne_in[ir])
  
    # Find fusion reaction rate
    # From the NRL plasma formulary
    source_in[ir] = ni2 * 3.68e-12 * Ti_kev^(-2.0/3.0) * exp(-19.94 * Ti_kev^(-1.0/3.0)) * 1.0e-6
  end
  
  # Then define the source in the computational domain
  for ir in 1:Nrad
    # Express Ti in kev
    Ti_kev = Ti[ir]/(1000.0*ev2joule)
   
    ni2 = (DTmix*ne[ir])*((1.0-DTmix)*ne[ir])
  
    # Find fusion reaction rate
    # From the NRL plasma formulary

    source[ir] = ni2 * 3.68e-12 * Ti_kev^(-2.0/3.0) * exp(-19.94 * Ti_kev^(-1.0/3.0)) * 1.0e-6
    source_local[ir] = source[ir]
    source[ir] = source[ir]*Vprime[ir]
  end
  
  for iv in 1:Nv
     idx = gindex(Nrad,iv)
     source[idx] = 0.0
     source_local[Nrad,iv] = 0.0
  end

  if zerosource
     source[:] = 0.0
     source_local[:,:] = 0.0
  end

end

function add2source(x)
  global source
  source = source+x
end

function add2source_element(idx,x)
  global source
  source[idx] = source[idx]+x
end


function zero_source_element(idx)
  global source
  source[idx] = 0.0
end

function set_source_element(idx,x)
  global source
  source[idx] = x
end

function gindex(ir,iv)
  return (ir-1)*Nv+iv
end

function init_source_for_analytic_test()
  global source

  source= zeros(Nrad*Nv)

  vhat = v/valpha
  rhat = rgrid/a

  for ir in 1:Nrad
    for iv in 1:Nv
      idx = gindex(ir,iv)

# Drr:
       source[idx] = 6.0*(2.0-rhat[ir])*v[iv]^2
# Dvv:
       source[idx] += 0.25*(2.0-rhat[ir]^2)*(6.0*vhat[iv]^2 - 2.0*vhat[iv]^3 - 4.0*vhat[iv]^4)*exp(-vhat[iv])*valpha^2
# Drv:
       source[idx] += brv*exp(-2.0*vhat[iv])*v[iv]^2 * ( (4.0*vhat[iv]/rhat[ir]) - 6.0*vhat[iv]*rhat[ir] )
# Dvr:
       source[idx] += bvr*4.0*exp(-2.0*vhat[iv])*rhat[ir]* ( valpha^2*vhat[iv] - v[iv]^2 - v[iv]^2*vhat[iv] )
       source[idx] *= D0*exp(-vhat[iv]^2)*Vprime[ir]/a^2

    end
  end
 
end

end
