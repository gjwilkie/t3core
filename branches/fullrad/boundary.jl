module boundary
using matrix: collop, gindex, find_local_sd, global_matrix, set_matrix_element, analytic_sd, broad_sd
using input
using sourcemod: add2source, source, set_source_element, source_local, add2source_element
using species: mass, Ti, ne
using grids: v, d3v, ddv, rgrid, ddr
using geometry: Vprime, surface_area
using matrix: global_matrix
using constants: valpha
using diffcoeff: Drr
using Dierckx
using PyPlot

export F0edge, calculate_boundary, fluxin

F0edge=Float64[]
fluxin=Float64[]


function maxwellian_f0(n,T)
   return n*(2.0*pi*T/m_trace)^(-1.5)*exp(-0.5*m_trace*v.^2/T)
end

function calculate_boundary()
  global F0edge, fluxin

  fluxin = zeros(Float64,Nv)
  if maxwellian_edge
    F0edge = maxwellian_f0(nedge,Ti[end]*Tashfac)
  else
    F0edge = broad_sd(Nrad,nedge,Ti[end]*Tashfac,false)
  end

  for iv in 1:Nv
    idx = gindex(Nrad,iv)
    for jdx in 1:Nv*Nrad
      set_matrix_element(idx,jdx,0.0)
    end
    set_matrix_element(idx,idx,1.0)
    set_source_element(idx,F0edge[iv]) 
  end

  # No incoming flux anymore

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
 
  Vp_func = Spline1D(rgrid_in,Vprime,k=spline_k)
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
