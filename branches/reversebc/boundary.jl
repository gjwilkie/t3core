module boundary
using matrix: collop, gindex, find_local_sd, global_matrix, set_matrix_element, analytic_sd
using input
using sourcemod: add2source, source, set_source_element, source_local, source_in, add2source_element
using species: mass, Ti, ne
using grids: v, d3v, ddv, rgrid, ddr
using geometry: Vprime, Vprime_global, surface_area, surface_area_global, grho_global
using matrix: global_matrix
using constants: valpha
using diffcoeff: Drr
using Dierckx

export F0edge, calculate_boundary, fluxout, fluxin

F0edge=Float64[]
fluxout=Float64[]
fluxin=Float64[]


function maxwellian_f0(n,T)
   return n*(2.0*pi*T/m_trace)^(-1.5)*exp(-0.5*m_trace*v.^2/T)
end

function calculate_boundary()
  global F0edge, fluxout

  fluxout = zeros(Float64,Nv)
  if maxwellian_edge
    F0edge = maxwellian_f0(nedge,Ti[1]*Tashfac)
  else
    F0edge = analytic_sd(1,nedge,Ti[1]*Tashfac,false)
  end

  for iv in 1:Nv
    idx = gindex(1,iv)
    for jdx in 1:Nv*Nrad
      set_matrix_element(idx,jdx,0.0)
    end
    set_matrix_element(idx,idx,1.0)
    set_source_element(idx,F0edge[iv]) 
  end

  # Get the rgrids between gs2/global right
#  integrand_func = Spline1D(rgrid_in,vec(source_in[:,1].*Vprime_global))

  integrated_source = zeros(Float64,Nv)
  for iv in 1:Nv
    function integrand(r)
  #    global integrand_func
      integrand_func = Spline1D(rgrid_in,vec(source_in[:,iv]).*Vprime_global,k=spline_k)
      return evaluate(integrand_func,r)
    end
    integrated_source[iv],err = quadgk(integrand,0.0,rgrid[end])
  end
 
  Vprime_func = Spline1D(rgrid_in,Vprime_global,k=spline_k)
  Vprime_h = evaluate(Vprime_func,rgrid[end] + 0.5*(rgrid[end]-rgrid[end-1]))

  tprimi = (Ti[end] - Ti[end-1])/(Ti[end]*(rgrid[end]-rgrid[end-1]))
  fprimi = (ne[end] - ne[end-1])/(ne[end]*(rgrid[end]-rgrid[end-1]))

  if ejection_mode
    # If this option is chosen, all alphas come in immediately from center
    fluxout = integrated_source/Vprime_h
    totalfluxout = dot(d3v,fluxout)
  else
    totalfluxout = dot(integrated_source,d3v)/Vprime_h
    # Otherwise, let what comes in be Maxwellian, such that one obtains total incoming particle flux when integrating over d3v

#    fluxin = totalfluxin*(m_trace/(2.0*pi*Ti[1]))^(1.5)*exp(-m_trace*v.^2/(2.0*Ti[1]))
#    fluxout = exp(-m_trace*v.^2/(2.0*Ti[end]*Tashfac)).*(fprimi + tprimi*( (0.5*m_trace*v.^2/(Ti[end]*Tashfac)) - 1.5) ).*vec(Drr[end,:])
    fluxin = exp(-m_trace*v.^2/(2.0*Ti[1]*Tashfac))

    normalize = dot(d3v,fluxout)
    fluxout = fluxout*totalfluxout/normalize
  end

  for iv in 1:Nv-1
    idx = gindex(Nrad,iv)
    add2source_element(idx, Vprime_h*fluxout[iv]*v[iv]^2/(rgrid[end]-rgrid[end-1]))
  end

end


function calculate_boundary_maxw()

  # Dirichlet BC at r=rmax: set density to be nedge (input), and temperature to be Ti @ rmax
  for ir in 1:2*Nrad
    set_matrix_element(Nrad,ir,0.0)
  end
  set_matrix_element(Nrad,Nrad,1.0)
  set_matrix_element(2*Nrad,2*Nrad,1.0)
  set_source_element(Nrad,nedge)
  set_source_element(2*Nrad,Ti[end])

  # Set the particle flux at r=rmin
  # No heat flux at r=rmin
  function integrand(r)
    return evaluate(integrand_func,r)
  end

  integrand_func = Spline1D(rgrid_in,vec(source_in[:]).*Vprime_global,k=spline_k)
  integrated_source,err = quadgk(integrand,0.0,rgrid_gs2[1])

  totalfluxin = integrated_source
  println("Total flux into domain = ",totalfluxin)

  area_func = Spline1D(rgrid_in,surface_area_global,k=spline_k)
  area_h = evaluate(area_func,rgrid[1] - 0.5*(rgrid[2]-rgrid[1]))

  add2source_element(1, area_h*totalfluxin/(rgrid[2]-rgrid[1]))

end

function calculate_boundary_for_analytic_test()
  global F0edge 
  F0edge = zeros(Float64,Nv)

  for iv in 1:Nv
    for ir in 1:Nrad
      idx = gindex(ir,iv)
      F0edge[iv]=exp(-v[iv]^2/valpha^2)
    end
    idx = gindex(Nrad,iv)
    set_source_element(idx,F0edge[iv])
  end 

#  rmin = (rgrid[1] - (rgrid[2] - rgrid[1]))/a
  rmin = (rgrid[1] - 0.5(rgrid[2] - rgrid[1]))/a
  vhat = v/valpha

  fluxin = Array(Float64,Nv)
  for iv in 1:Nv
    # Flux = - Drr * df/dr - Drv * df/dv
    fluxin[iv] = (D0*(3.0-rmin))*2.0*(rmin/a)*exp(-vhat[iv]^2)
#    if iv > 1
      fluxin[iv] = fluxin[iv] + (D0*brv*exp(-2.0*vhat[iv])*valpha/a) * (2.0-rmin^2) * (2*vhat[iv]/valpha) * exp(-vhat[iv]^2)
#    end
  end
 
  Vp_func = Spline1D(rgrid_in,Vprime_global,k=spline_k)
  Vp_h = evaluate(Vp_func,rgrid[1] - 0.5*(rgrid[2]-rgrid[1]))

#  for iv in 2:Nv-1
  for iv in 1:Nv-1
    idx = gindex(1,iv)
    add2source_element(idx, Vp_h*fluxin[iv]*v[iv]^2/(rgrid[2]-rgrid[1]))
  end

#  for ir in 1:Nrad
#    idx = gindex(ir,1)
#    set_source_element(idx,0.0)
#  end

end
 
end
