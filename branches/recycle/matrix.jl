module matrix
using input: Nrad, Nv, rgrid_gs2, tracespecs, Tashfac, nedge, deltat, m_trace, Z_trace, DTmix, semianalytic_on, rgrid_in, spline_k, recycle
using collisions: lnLambda
using constants
using diffcoeff: Drr,Drv,Dvr,Dvv, Dnn, DnT, DTn, DTT, pflux0, hflux0
using sourcemod: source, zero_source_element, source_local, reaction_rate, add2source_element, set_source_element
using chebyshev
using geometry: Vprime, surface_area, surface_area_global, Vprime_global
using species: Te, ne, Ti, vcrit
using grids: init_vgrid, v, d3v, ddv, ddr, init_rgrid, rgrid
using Base.Test
using Dierckx

export build_matrix, init_collop, solve_steadystate, f0, gindex, vindex, rindex, C_global, ddr_global, ddv_global, Dv_global, collop, global_matrix, Vprime_global, v_global, nupar, taus, set_matrix_element, collop_ion, zero_collop, nu_s_v3, nu_par_v3, collop_el, build_matrix_maxw

taus = Array(Float64,1)
collop = Array(Float64,3)
collop_ion = Array(Float64,3)
collop_el = Array(Float64,3)
f0 = Array(Float64,1)
global_matrix= Array(Float64,2)
nupar= Array(Float64,2)
nu_s_v3 = Array(Float64,2)
nu_par_v3 = Array(Float64,2)

# Inializes global matrix

function build_matrix()
  global collop, global_matrix, nu_s_v3, nu_par_v3


  # Matrix equation solves for g = F0(r,v) - Fedge(v)
  # subject to the boundary conditions dg/dr=0 at r=0 and g=0 and r=a
 
  global_matrix=zeros(Nrad*Nv,Nrad*Nv)

  # Prepare arrays for easier syntax
  delta_r = rgrid[2] - rgrid[1]
  delta_v = v[2] - v[1]

  v_jph = 0.5*collect(v[1:end-1]+v[2:end])
  v_jmh = copy(v_jph)
  push!(v_jph,v[end]+0.5*(v[end]-v[end-1]))
  unshift!(v_jmh,v[1] - 0.5(v[2]-v[1]))

  r_iph = 0.5*collect(rgrid[1:end-1]+rgrid[2:end])
  r_imh = copy(r_iph)
  push!(r_iph,rgrid[end]+0.5*(rgrid[end]-rgrid[end-1]))
  unshift!(r_imh,rgrid[1] - 0.5(rgrid[2]-rgrid[1]))

  Vprime_func = Spline1D(rgrid_in,Vprime_global,k=spline_k)
  Vprime_iph = evaluate(Vprime_func,r_iph)
  Vprime_imh = evaluate(Vprime_func,r_imh)

  Drr_iph = zeros(Nrad,Nv)
  Drr_imh = zeros(Nrad,Nv)
  Drv_iph = zeros(Nrad,Nv)
  Drv_imh = zeros(Nrad,Nv)
  for jv in 1:Nv
    Drr_func = Spline1D(rgrid,vec(Drr[:,jv]),k=spline_k)
    Drr_iph[:,jv] = evaluate(Drr_func,r_iph)
    Drr_imh[:,jv] = evaluate(Drr_func,r_imh)
    Drv_func = Spline1D(rgrid,vec(Drv[:,jv]),k=spline_k)
    Drv_iph[:,jv] = evaluate(Drv_func,r_iph)
    Drv_imh[:,jv] = evaluate(Drv_func,r_imh)
  end
  nus_jph = zeros(Nrad,Nv)
  nus_jmh = zeros(Nrad,Nv)
  nupar_jph = zeros(Nrad,Nv)
  nupar_jmh = zeros(Nrad,Nv)
  Dvr_jph = zeros(Nrad,Nv)
  Dvr_jmh = zeros(Nrad,Nv)
  Dvv_jph = zeros(Nrad,Nv)
  Dvv_jmh = zeros(Nrad,Nv)
  for ir in 1:Nrad
    nus_func = Spline1D(v,vec(nu_s_v3[ir,:]),k=1)
    nus_jph[ir,:] = evaluate(nus_func,v_jph)'
    nus_jmh[ir,:] = evaluate(nus_func,v_jmh)'
    nupar_func = Spline1D(v,vec(nu_par_v3[ir,:]),k=1)
    nupar_jph[ir,:] = evaluate(nupar_func,v_jph)'
    nupar_jmh[ir,:] = evaluate(nupar_func,v_jmh)'
    Dvr_func = Spline1D(v,vec(Dvr[ir,:]),k=1)
    Dvr_jph[ir,:] = evaluate(Dvr_func,v_jph)'
    Dvr_jmh[ir,:] = evaluate(Dvr_func,v_jmh)'
    Dvv_func = Spline1D(v,vec(Dvv[ir,:]),k=1)
    Dvv_jph[ir,:] = evaluate(Dvv_func,v_jph)'
    Dvv_jmh[ir,:] = evaluate(Dvv_func,v_jmh)'
  end
  
  nupar_jmh[:,1] = 0.0
  nus_jmh[:,1] = 0.0
  
#  for iv in 1:Nv
#    println("nus = ",nus_jph[:,iv]/v[iv]^3)
#    println("nupar = ",nus_jph[:,iv]/v[iv]^3)
#  end
 
  for idx in 1:Nrad*Nv
    ir = rindex(idx)
    jv= vindex(idx)
  
    # Build stencil matrices
    flux_rr_iph = zeros(Nrad*Nv)
    flux_rr_imh = zeros(Nrad*Nv)
    flux_rv_iph = zeros(Nrad*Nv)
    flux_rv_imh = zeros(Nrad*Nv)
    flux_vr_jph = zeros(Nrad*Nv)
    flux_vr_jmh = zeros(Nrad*Nv)
    flux_vv_jph = zeros(Nrad*Nv)
    flux_vv_jmh = zeros(Nrad*Nv)
    nus_term    = zeros(Nrad*Nv)
    nupar_term  = zeros(Nrad*Nv)


    # Build the various terms in the steady-state equation. 
  
    # Internal points
    if (2 <= ir <= Nrad-1) && (2 <= jv <= Nv-1)
      flux_rr_iph[gindex(ir+1,jv)] += (Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_iph[gindex(ir,jv)] += -(Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir,jv)] += (Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir-1,jv)] += -(Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r

      flux_rv_iph[gindex(ir+1,jv+1)] += (Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir+1,jv-1)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir,jv+1)] += (Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir,jv-1)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir,jv-1)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir-1,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir-1,jv-1)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)

      flux_vv_jph[gindex(ir,jv+1)] += (Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv)] += (Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv-1)] += -(Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v

      flux_vr_jph[gindex(ir+1,jv+1)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jph[gindex(ir-1,jv+1)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jph[gindex(ir+1,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jph[gindex(ir-1,jv)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jmh[gindex(ir+1,jv)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
      flux_vr_jmh[gindex(ir-1,jv)] += -(Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
      flux_vr_jmh[gindex(ir+1,jv-1)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
      flux_vr_jmh[gindex(ir-1,jv-1)] += -(Dvr_jmh[ir,jv] *Vprime[ir]* v_jmh[jv]^2) * (0.25/delta_r)

      nupar_term[gindex(ir,jv+1)] += 0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += - 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv-1)] += 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv+1)] += 0.5*nus_jph[ir,jv]*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv]-nus_jmh[ir,jv])*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv-1)] += -0.5*nus_jmh[ir,jv]*Vprime[ir] /(delta_v)
    # f = 0 at v=vmax, so remove all coefficients that would multiply such a point
    elseif (2 <= ir <= Nrad-1) && (jv == Nv) 
      flux_rr_iph[gindex(ir+1,jv)] += (Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_iph[gindex(ir,jv)] += -(Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir,jv)] += (Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir-1,jv)] += -(Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r

      flux_rv_iph[gindex(ir+1,jv-1)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir,jv-1)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir,jv-1)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir-1,jv-1)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)

      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv)] += (Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv-1)] += -(Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v

      flux_vr_jph[gindex(ir+1,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jph[gindex(ir-1,jv)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jmh[gindex(ir+1,jv)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
      flux_vr_jmh[gindex(ir-1,jv)] += -(Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
      flux_vr_jmh[gindex(ir+1,jv-1)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
      flux_vr_jmh[gindex(ir-1,jv-1)] += -(Dvr_jmh[ir,jv] *Vprime[ir]* v_jmh[jv]^2) * (0.25/delta_r)

      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += - 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv-1)] += 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv]-nus_jmh[ir,jv])*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv-1)] += -0.5*nus_jmh[ir,jv]*Vprime[ir] /(delta_v)

    # df/dv = 0 at v=0, in the sense that "v[0]" = v[2], all centered differences vanish
    elseif (2 <= ir <= Nrad-1) && (jv == 1)
      flux_rr_iph[gindex(ir+1,jv)] += (Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_iph[gindex(ir,jv)] += -(Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir,jv)] += (Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir-1,jv)] += -(Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r

      flux_rv_iph[gindex(ir+1,jv+1)] += (Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_iph[gindex(ir+1,jv)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_iph[gindex(ir,jv+1)] += (Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_iph[gindex(ir,jv)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_imh[gindex(ir,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_imh[gindex(ir,jv)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_imh[gindex(ir-1,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_imh[gindex(ir-1,jv)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)

      flux_vv_jph[gindex(ir,jv+1)] += (Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
#      flux_vv_jmh[gindex(ir,jv)] += (Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v
#      flux_vv_jmh[gindex(ir,jv-1)] += -(Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v

      flux_vr_jph[gindex(ir+1,jv+1)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jph[gindex(ir-1,jv+1)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jph[gindex(ir+1,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
      flux_vr_jph[gindex(ir-1,jv)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir+1,jv)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir-1,jv)] += -(Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir+1,jv-1)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir-1,jv-1)] += -(Dvr_jmh[ir,jv] *Vprime[ir]* v_jmh[jv]^2) * (0.25/delta_r)


      nupar_term[gindex(ir,jv+1)] += 0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv+1)] += 0.5*nus_jph[ir,jv]*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv])*Vprime[ir] /(delta_v)


    # Flux specified at r=rmin
    elseif ( ir == 1) && (2 <= jv <= Nv-1)
      flux_rr_iph[gindex(ir+1,jv)] += (Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_iph[gindex(ir,jv)] += -(Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r

      flux_rv_iph[gindex(ir+1,jv+1)] += (Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir+1,jv-1)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir,jv+1)] += (Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir,jv-1)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)

      flux_vv_jph[gindex(ir,jv+1)] += (Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv)] += (Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv-1)] += -(Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v

      flux_vr_jph[gindex(ir+1,jv+1)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir,jv+1)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir+1,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir,jv)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir+1,jv)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir,jv)] += -(Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir+1,jv-1)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir,jv-1)] += -(Dvr_jmh[ir,jv] *Vprime[ir]* v_jmh[jv]^2) * (0.5/delta_r)

      nupar_term[gindex(ir,jv+1)] += 0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += - 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv-1)] += 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv+1)] += 0.5*nus_jph[ir,jv]*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv]-nus_jmh[ir,jv])*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv-1)] += -0.5*nus_jmh[ir,jv]*Vprime[ir] /(delta_v)


    # Recycling at r=rmax
    elseif ( ir == Nrad) && (2 <= jv <= Nv-1)
      flux_rr_imh[gindex(ir,jv)] += (Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir-1,jv)] += -(Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r

      flux_rr_iph[gindex(ir,jv)] += (1.0-recycle)*flux_rr_imh[gindex(ir,jv)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rr_iph[gindex(ir-1,jv)] += (1.0-recycle)*flux_rr_imh[gindex(ir-1,jv)]*Vprime_iph[ir]/Vprime_imh[ir]

      flux_rv_imh[gindex(ir,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir,jv-1)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir-1,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir-1,jv-1)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir,jv+1)] += (1.0-recycle)*flux_rv_imh[gindex(ir,jv+1)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rv_iph[gindex(ir,jv-1)] += (1.0-recycle)*flux_rv_imh[gindex(ir,jv-1)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rv_iph[gindex(ir-1,jv+1)] += (1.0-recycle)*flux_rv_imh[gindex(ir-1,jv+1)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rv_iph[gindex(ir-1,jv-1)] += (1.0-recycle)*flux_rv_imh[gindex(ir-1,jv-1)]*Vprime_iph[ir]/Vprime_imh[ir]

      flux_vv_jph[gindex(ir,jv+1)] += (Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv)] += (Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv-1)] += -(Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v

      flux_vr_jph[gindex(ir,jv+1)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir-1,jv+1)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir-1,jv)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir,jv)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir-1,jv)] += -(Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir,jv-1)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir-1,jv-1)] += -(Dvr_jmh[ir,jv] *Vprime[ir]* v_jmh[jv]^2) * (0.5/delta_r)

      nupar_term[gindex(ir,jv+1)] += 0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += - 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv-1)] += 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv+1)] += 0.5*nus_jph[ir,jv]*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv]-nus_jmh[ir,jv])*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv-1)] += -0.5*nus_jmh[ir,jv]*Vprime[ir] /(delta_v)
 
    elseif ( ir == Nrad) && (jv == 1)
      flux_rr_imh[gindex(ir,jv)] += (Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir-1,jv)] += -(Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
      flux_rr_iph[gindex(ir,jv)] += (1.0-recycle)*flux_rr_imh[gindex(ir,jv)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rr_iph[gindex(ir-1,jv)] += (1.0-recycle)*flux_rr_imh[gindex(ir-1,jv)]*Vprime_iph[ir]/Vprime_imh[ir]

      flux_rv_imh[gindex(ir,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_imh[gindex(ir,jv)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_imh[gindex(ir-1,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_imh[gindex(ir-1,jv)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_iph[gindex(ir,jv+1)] += (1.0-recycle)*flux_rv_imh[gindex(ir,jv+1)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rv_iph[gindex(ir,jv)] += (1.0-recycle)*flux_rv_imh[gindex(ir,jv)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rv_iph[gindex(ir-1,jv+1)] += (1.0-recycle)*flux_rv_imh[gindex(ir-1,jv+1)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rv_iph[gindex(ir-1,jv)] += (1.0-recycle)*flux_rv_imh[gindex(ir-1,jv)]*Vprime_iph[ir]/Vprime_imh[ir]

      flux_vv_jph[gindex(ir,jv+1)] += (Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v

      flux_vr_jph[gindex(ir,jv+1)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir-1,jv+1)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir-1,jv)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)


      nupar_term[gindex(ir,jv+1)] += 0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv+1)] += 0.5*nus_jph[ir,jv]*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv])*Vprime[ir] /(delta_v)


    elseif ( ir == Nrad) && (jv == Nv)
      flux_rr_imh[gindex(ir,jv)] += (Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
      flux_rr_imh[gindex(ir-1,jv)] += -(Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
      flux_rr_iph[gindex(ir,jv)] += (1.0-recycle)*flux_rr_imh[gindex(ir,jv)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rr_iph[gindex(ir-1,jv)] += (1.0-recycle)*flux_rr_imh[gindex(ir-1,jv)]*Vprime_iph[ir]/Vprime_imh[ir]

      flux_rv_imh[gindex(ir,jv-1)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_imh[gindex(ir-1,jv-1)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir,jv-1)] += (1.0-recycle)*flux_rv_imh[gindex(ir,jv-1)]*Vprime_iph[ir]/Vprime_imh[ir]
      flux_rv_iph[gindex(ir-1,jv-1)] += (1.0-recycle)*flux_rv_imh[gindex(ir-1,jv-1)]*Vprime_iph[ir]/Vprime_imh[ir]

      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv)] += (Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv-1)] += -(Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v

      flux_vr_jph[gindex(ir,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir-1,jv)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir,jv)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir-1,jv)] += -(Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir,jv-1)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.5/delta_r)
      flux_vr_jmh[gindex(ir-1,jv-1)] += -(Dvr_jmh[ir,jv] *Vprime[ir]* v_jmh[jv]^2) * (0.5/delta_r)

      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += - 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv-1)] += 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv]-nus_jmh[ir,jv])*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv-1)] += -0.5*nus_jmh[ir,jv]*Vprime[ir] /(delta_v)

    # r=0 Corner cases:
    elseif (ir == 1) && (jv == 1)
      flux_rr_iph[gindex(ir+1,jv)] += (Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_iph[gindex(ir,jv)] += -(Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
#      flux_rr_imh[gindex(ir,jv)] += (Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r
#      flux_rr_imh[gindex(ir-1,jv)] += -(Drr_imh[ir,jv]*Vprime_imh[ir]*v[jv]^2)/delta_r

      flux_rv_iph[gindex(ir+1,jv+1)] += (Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_iph[gindex(ir+1,jv)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_iph[gindex(ir,jv+1)] += (Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.5/delta_v)
      flux_rv_iph[gindex(ir,jv)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.5/delta_v)
#      flux_rv_imh[gindex(ir,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
#      flux_rv_imh[gindex(ir,jv)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
#      flux_rv_imh[gindex(ir-1,jv+1)] += (Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)
#      flux_rv_imh[gindex(ir-1,jv)] += -(Drv_imh[ir,jv]*Vprime_imh[ir] * v[jv]^2) * (0.5/delta_v)

      flux_vv_jph[gindex(ir,jv+1)] += (Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
#      flux_vv_jmh[gindex(ir,jv)] += (Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v
#      flux_vv_jmh[gindex(ir,jv-1)] += -(Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v

      flux_vr_jph[gindex(ir+1,jv+1)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir,jv+1)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir+1,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
      flux_vr_jph[gindex(ir,jv)] += -(Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.5/delta_r)
#      flux_vr_jmh[gindex(ir+1,jv)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir-1,jv)] += -(Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir+1,jv-1)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir-1,jv-1)] += -(Dvr_jmh[ir,jv] *Vprime[ir]* v_jmh[jv]^2) * (0.25/delta_r)


      nupar_term[gindex(ir,jv+1)] += 0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv+1)] += 0.5*nus_jph[ir,jv]*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv])*Vprime[ir] /(delta_v)

    elseif ( ir == 1) && (jv == Nv)
      flux_rr_iph[gindex(ir+1,jv)] += (Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r
      flux_rr_iph[gindex(ir,jv)] += -(Drr_iph[ir,jv]*Vprime_iph[ir]*v[jv]^2)/delta_r

      flux_rv_iph[gindex(ir+1,jv-1)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)
      flux_rv_iph[gindex(ir,jv-1)] += -(Drv_iph[ir,jv]*Vprime_iph[ir] * v[jv]^2) * (0.25/delta_v)

      flux_vv_jph[gindex(ir,jv)] += -(Dvv_jph[ir,jv]*Vprime[ir]*v_jph[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv)] += (Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v
      flux_vv_jmh[gindex(ir,jv-1)] += -(Dvv_jmh[ir,jv]*Vprime[ir]*v_jmh[jv]^2)/delta_v

#      flux_vr_jph[gindex(ir+1,jv)] += (Dvr_jph[ir,jv]*Vprime[ir] * v_jph[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir+1,jv)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)
#      flux_vr_jmh[gindex(ir+1,jv-1)] += (Dvr_jmh[ir,jv]*Vprime[ir] * v_jmh[jv]^2) * (0.25/delta_r)

      nupar_term[gindex(ir,jv)] += -0.5*nupar_jph[ir,jv]*v_jph[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv)] += - 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nupar_term[gindex(ir,jv-1)] += 0.5*nupar_jmh[ir,jv]*v_jmh[jv]*Vprime[ir]/(delta_v^2)
      nus_term[gindex(ir,jv)] += 0.5*(nus_jph[ir,jv]-nus_jmh[ir,jv])*Vprime[ir] /(delta_v)
      nus_term[gindex(ir,jv-1)] += -0.5*nus_jmh[ir,jv]*Vprime[ir] /(delta_v)

    else
      error("Failed to catch idx as a boundary condition.")
    end

    global_matrix[idx,:] += - ( flux_rr_iph - flux_rr_imh + flux_rv_iph - flux_rv_imh )' / delta_r 
    global_matrix[idx,:] += - ( flux_vr_jph - flux_vr_jmh + flux_vv_jph - flux_vv_jmh )' / delta_v
    if !semianalytic_on
      global_matrix[idx,:] += - ( nus_term + nupar_term )'
    end

  end

end 

function init_collop()
  global collop, nupar, taus, collop_ion, nu_s_v3, nu_par_v3, collop_el
  
  # Collision operator is a linear Nv-by-Nv operator acting on a distribution function
  collop = zeros(Nv,Nv,Nrad)
  collop_ion = zeros(Nv,Nv,Nrad)
  collop_el = zeros(Nv,Nv,Nrad)
  taus=zeros(Nrad)
  nupar = zeros(Nrad,Nv)
  nu_s_v3 = zeros(Nrad,Nv)
  nu_par_v3 = zeros(Nrad,Nv)

  # Find the collision operator against the static background at each radius
  for ir in 1:Nrad
    # Coulomb logarithm
    # From pg. 34 of NRL formulary
    # Should probably be different for each species, but this is more complicated
    # than expected and won't make that much difference. For now, just using the value
    # for electrons.
    T_coll = [Ti[ir], Ti[ir], Te[ir]]
    n_coll = [DTmix*ne[ir], (1.0-DTmix)*ne[ir], ne[ir]]
    m_coll = [2.0*mp, 3.0*mp, me]
    Z_coll = [1.0,1.0,-1.0]
 
    # Collision operator assumes a mixture of D-T and electrons
    for is in 1:3
      # Coulomb logarithm assumes T_alpha = 10*T_i and n_alpha = 0.01*ne. Better ideas?
      logLambda = lnLambda(m1=m_trace, m2=m_coll[is], Z1=Z_trace, Z2=Z_coll[is], T1=10.0*Ti[ir], T2=T_coll[is], n1=0.01*ne[ir], n2=n_coll[is])

      vts = sqrt(2.0*T_coll[is]/m_coll[is])

      G = erf(v/vts) - (2.0/sqrt(pi))*(v/vts).*exp(-(v/vts).^2)
      G = G./(2.0*(v/vts).^2)

      # Construct (nuhat*vta^3) for this species pair
      # Eq. (3.48) of Hellander and Sigmar, times vta^3
      nuhat_vta3 = n_coll[is]*Z_coll[is]^2*Z_trace^2*el^4*logLambda
      nuhat_vta3 = nuhat_vta3 / (4.0*pi*(ep0*m_trace)^2)

      nu_s_v3[ir,:] = nu_s_v3[ir,:] + nuhat_vta3*m_trace*(G.*v.*v)[:]'/T_coll[is]

      nu_par_v3[ir,:] = nu_par_v3[ir,:] + 2.0*nuhat_vta3*G[:]'

      nupar[ir,:] = nupar[ir,:] + nu_par_v3[ir,:]./(v[:].^3)'

      if (is == 2)
        tempMatrix= diagm(vec(nu_s_v3[ir,:])) + diagm(0.5*v)*diagm(vec(nu_par_v3[ir,:]))*ddv
        collop_ion[:,:,ir] = diagm(v.^(-2))*(ddv*tempMatrix)
      end
      if (is == 3)
        tempMatrix= diagm(vec(nuhat_vta3*m_trace*(G.*v.*v)'/T_coll[is])) + diagm(v)*diagm(vec(nuhat_vta3*G'))*ddv
        collop_el[:,:,ir] = diagm(v.^(-2))*(ddv*tempMatrix)
      end 
        
    end

    tempMatrix= diagm(vec(nu_s_v3[ir,:])) + diagm(0.5*v)*diagm(vec(nu_par_v3[ir,:]))*ddv

    collop[:,:,ir] = diagm(v.^(-2))*(ddv*tempMatrix)

    # Regularlize so that df/dv=0 at v=0
#    collop[1,:,ir] = ddv[1,:] 
#    collop_ion[1,:,ir] = ddv[1,:] 
#    collop_el[1,:,ir] = ddv[1,:] 

    logLambda = lnLambda(m1=m_trace, m2=me, Z1=Z_trace, Z2=-1, T1=10.0*Ti[ir], T2=Te[ir], n1=0.01*ne[ir], n2=ne[ir])
    taus[ir] = (3.0/(16.0*sqrt(pi)))*me*m_trace*sqrt(2.0*Te[ir]/me)^3*(4.0*pi*ep0)^2/(Z_trace^2*el^4*ne[ir]*logLambda)
  end

end

function build_matrix_maxw()
  global global_matrix

  global_matrix = zeros(Float64,(2*Nrad,2*Nrad))
 
  delta_r = rgrid[2] - rgrid[1]

  r_iph = 0.5*(rgrid[1:end-1]+rgrid[2:end])
  r_imh = copy(r_iph)
  push!(r_iph,rgrid[end]+0.5*(rgrid[end]-rgrid[end-1]))
  unshift!(r_imh,rgrid[1] - 0.5(rgrid[2]-rgrid[1]))

  Vprime_func = Spline1D(rgrid,surface_area)
  Vprime_iph = evaluate(Vprime_func,r_iph)
  Vprime_imh = evaluate(Vprime_func,r_imh)

  Dnn_iph = zeros(Nrad)
  Dnn_imh = zeros(Nrad)
  DnT_iph = zeros(Nrad)
  DnT_imh = zeros(Nrad)
  DTn_iph = zeros(Nrad)
  DTn_imh = zeros(Nrad)
  DTT_iph = zeros(Nrad)
  DTT_imh = zeros(Nrad)
  pflux0_imh = zeros(Nrad)
  pflux0_iph = zeros(Nrad)
  hflux0_imh = zeros(Nrad)
  hflux0_iph = zeros(Nrad)

  Dnn_func = Spline1D(rgrid,vec(Dnn[:]))
  Dnn_iph[:] = evaluate(Dnn_func,r_iph)
  Dnn_imh[:] = evaluate(Dnn_func,r_imh)
  DnT_func = Spline1D(rgrid,vec(DnT[:]))
  DnT_iph[:] = evaluate(DnT_func,r_iph)
  DnT_imh[:] = evaluate(DnT_func,r_imh)
  DTn_func = Spline1D(rgrid,vec(DTn[:]))
  DTn_iph[:] = evaluate(DTn_func,r_iph)
  DTn_imh[:] = evaluate(DTn_func,r_imh)
  DTT_func = Spline1D(rgrid,vec(DTT[:]))
  DTT_iph[:] = evaluate(DTT_func,r_iph)
  DTT_imh[:] = evaluate(DTT_func,r_imh)
  pflux0_func = Spline1D(rgrid,vec(pflux0[:]))
  pflux0_iph[:] = evaluate(pflux0_func,r_iph)
  pflux0_imh[:] = evaluate(pflux0_func,r_imh)
  hflux0_func = Spline1D(rgrid,vec(hflux0[:]))
  hflux0_iph[:] = evaluate(hflux0_func,r_iph)
  hflux0_imh[:] = evaluate(hflux0_func,r_imh)

  for ir in 2:Nrad-1
    # Particle flux at zero gradient
    add2source_element(ir,-(Vprime_iph[ir]*pflux0_iph[ir]-Vprime_imh[ir]*pflux0_imh[ir])/delta_r)

    # Particle flux due to density gradient
    global_matrix[ir,ir-1] += -Vprime_imh[ir]*Dnn_imh[ir] /delta_r^2
    global_matrix[ir,ir] += (Vprime_imh[ir]*Dnn_imh[ir]+Vprime_iph[ir]*Dnn_iph[ir]) /delta_r^2
    global_matrix[ir,ir+1] += -Vprime_iph[ir]*Dnn_iph[ir] /delta_r^2

    # Particle flux due to temp gradient
    global_matrix[ir,Nrad+ir-1] += -Vprime_imh[ir]*DnT_imh[ir] /delta_r^2
    global_matrix[ir,Nrad+ir] += (Vprime_imh[ir]*DnT_imh[ir]+Vprime_iph[ir]*DnT_iph[ir]) /delta_r^2
    global_matrix[ir,Nrad+ir+1] += -Vprime_iph[ir]*DnT_iph[ir] /delta_r^2

    # Heat flux at zero gradient
    add2source_element(Nrad+ir,-(Vprime_iph[ir]*hflux0_iph[ir]-Vprime_imh[ir]*hflux0_imh[ir])/delta_r)

    # Heat flux due to density gradient
    global_matrix[Nrad+ir,ir-1] += -Vprime_imh[ir]*DTn_imh[ir] /delta_r^2
    global_matrix[Nrad+ir,ir] += (Vprime_imh[ir]*DTn_imh[ir]+Vprime_iph[ir]*DTn_iph[ir]) /delta_r^2
    global_matrix[Nrad+ir,ir+1] += -Vprime_iph[ir]*DTn_iph[ir] /delta_r^2

    # Heat flux due to temperature gradient
    global_matrix[Nrad+ir,Nrad+ir-1] += -Vprime_imh[ir]*DTT_imh[ir] /delta_r^2
    global_matrix[Nrad+ir,Nrad+ir] += (Vprime_imh[ir]*DTT_imh[ir]+Vprime_iph[ir]*DTT_iph[ir]) /delta_r^2
    global_matrix[Nrad+ir,Nrad+ir+1] += -Vprime_iph[ir]*DTT_iph[ir] /delta_r^2
  end

  # Boundary condition at rmin. Flux at minus location gets put to source instead
  # Particle flux at zero gradient
  add2source_element(1,-(Vprime_iph[1]*pflux0_iph[1])/delta_r)

  # Particle flux due to density gradient
  global_matrix[1,1] += Vprime_iph[1]*Dnn_iph[1] /delta_r^2
  global_matrix[1,2] += -Vprime_iph[1]*Dnn_iph[1] /delta_r^2

  # Particle flux due to temp gradient
  global_matrix[1,Nrad+1] += Vprime_iph[1]*DnT_iph[1] /delta_r^2
  global_matrix[1,Nrad+2] += -Vprime_iph[1]*DnT_iph[1] /delta_r^2

  # Heat flux at zero gradient
  add2source_element(Nrad+1,-(Vprime_iph[1]*hflux0_iph[1])/delta_r)

  # Heat flux due to density gradient
  global_matrix[Nrad+1,1] += Vprime_iph[1]*DTn_iph[1] /delta_r^2
  global_matrix[Nrad+1,2] += -Vprime_iph[1]*DTn_iph[1] /delta_r^2

  # Heat flux due to temperature gradient
  global_matrix[Nrad+1,Nrad+1] += Vprime_iph[1]*DTT_iph[1] /delta_r^2
  global_matrix[Nrad+1,Nrad+2] += -Vprime_iph[1]*DTT_iph[1] /delta_r^2

  # Boundary condition at r=rmax is easy:
  global_matrix[Nrad,:] = 0.0
  global_matrix[Nrad,Nrad] = 1.0
  set_source_element(Nrad, nedge)
  global_matrix[2*Nrad,:] = 0.0
  global_matrix[2*Nrad,2*Nrad] = 1.0
  set_source_element(2*Nrad,Ti[end])

end

function solve_steadystate()
  global f0, Vprime_global, v_global
  # Equation was multipled by V'(rho) to avoid division by zero, so here we make sure the source term reflects that
 
  f0 = global_matrix\source
#  f0 = pinv(global_matrix)*source
end

function advance_timestep(it)
  global f0, global_matrix 
  if it ==1
    f0 = zeros(Nrad*Nv)
  end
  f0 = (eye(Nrad*Nv)+deltat*global_matrix)\(source*deltat+f0)
end

function find_local_sd(ir,nlocal)
  global collop
  F0local = zeros(Float64,Nv+1)
  matrix = zeros(Float64,Nv+1,Nv+1)
  RHS = zeros(Nv+1)

  for iv in 1:Nv
    for jv in 1:Nv
      matrix[iv,jv] = -collop[iv,jv,ir]
    end
    RHS[iv] = source_local[ir,iv]
  end
  RHS[Nv+1] = nlocal

  matrix[1,:] = 0.0
  matrix[1,1:Nv] = ddv[1,:]
  RHS[1] = 0.0

  sinkdist = exp(-m_trace*v.^2/(2.0*Ti[ir]))
  matrix[2:Nv,Nv+1] = sinkdist[2:Nv]

  matrix[Nv+1,1:Nv] = d3v[:]

#  matrix[Nv+1,:] = matrix[Nv+1,:]/nedge
#  RHS[Nv+1] = RHS[Nv+1]/nedge

  F0local = matrix\RHS

  return F0local[1:Nv]

end

function analytic_sd(ir,n_in,Tash_in,use_as_nash)
  global taus

  fsd = (reaction_rate[ir]*taus[ir]*0.25/pi)./(vcrit[ir]^3 + v.^3)
  idx = indmin(abs(v-valpha))
  fsd[idx:end] = 0.0

  nsd = dot(d3v,fsd)

  if use_as_nash
    # Use the input density as the ash density
    nash_use = n_in
  else
    # Use the input density as the total alpha density
    nash_use = n_in - nsd
  end

  ashdist = exp(-m_trace*v.^2/(2.0*Tash_in))
  ashdist = ashdist*nash_use*(pi*2.0*Tash_in/m_trace)^(-1.5)

  return fsd + ashdist
end

function set_matrix_element(idx,jdx,x)
  global global_matrix
  global_matrix[idx,jdx]  = x
end

function zero_collop()
  global collop, collop_ion, nu_s_v3, nu_par_v3, nupar
  collop = zeros(Nv,Nv,Nrad)
  collop_ion = zeros(Nv,Nv,Nrad)
  nu_s_v3 = zeros(Nrad,Nv)
  nu_par_v3 = zeros(Nrad,Nv)
  nupar = zeros(Nrad,Nv)
  for ir in 1:Nrad
    for iv in 1:Nv
      nu_par_v3[ir,iv] = 0.0* v[iv]^3/valpha^3
      nu_s_v3[ir,iv] = nu_par_v3[ir,iv] * v[iv]^2/valpha^2
    end
  end
end

# Finds the global matrix index given the velocity and radial indices
function gindex(ir,iv)
  return (ir-1)*Nv+iv
end

# Finds the radial index given a global index
function rindex(idx)
  return div(idx-1,Nv)+1
end

# Finds the velocity index given a global index
function vindex(idx)
  return mod(idx-1,Nv)+1
end

end

