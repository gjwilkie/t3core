module postproc
#using Winston
using grids: v, d3v, rgrid, ddv,ddr
using input: Nv,Nrad,Nt, deltat, diffmodel, a, tracespecs, m_trace, ir_sample, rgrid_gs2, rhostar, mref, rmaj, vmax, DTmix, Z_trace, Te_in, Ti_in, ne_in, rgrid_in, nedge, turbfac, Tashfac,vt_temp_fac, plot_output, adj_postproc_temp
using matrix: f0, gindex, nupar, find_local_sd, analytic_sd, collop_ion, collop, nu_par_v3, nu_s_v3, collop_el, taus, broad_sd, init_collop
using diffcoeff: Drv, Drr, Dvr, Dvv, chii,phi2,hflux_tot
using geometry: surface_area_global, Vprime, grad_rho
using species: mass, Ti, ne, nref, Tref, Te, vcrit
using constants: valpha, Ealpha, ep0, mp, me, el
using boundary: F0edge, fluxin
using sourcemod: source_local, reaction_rate, source_in
using collisions: G, lnLambda
using PyPlot
using CurveFit
using Cubature
using Dierckx

export plot_steadystate, f0sd, transient_diagnostics

f0save = Array(Float64,3)
ir_set = Int64[]
basedir = ASCIIString[]
dir = ASCIIString[]
Tash = Float64[]
nash = Float64[]

function init_postproc(runname)
  global basedir, dir
  basedir = copy(runname)
  if deltat < 0.0
    dir = copy(basedir)
  end

end

function transient_diagnostics(f0in,it)
   global dir
   dir = basedir*"/t"*string(it)
   try run(`mkdir $dir`)
   catch
     println("Directory "*dir*"already exists")
   end
 
   plot_steadystate(f0in)

#   println(f0in[(ir_sample-1)*Nv+1:ir_sample*Nv])
   semilogy(v/valpha,abs(vec(f0in[(ir_sample-1)*Nv+1:ir_sample*Nv])))

end

function plot_steadystate(f0in)
  global dir, ir_set, Tash,nash

  f0alpha = zeros(Float64,(Nrad,Nv))
  for ir in 1:Nrad
     f0alpha[ir,:] = vec(f0in[(ir-1)*Nv+1:ir*Nv])
  end

  ashfit = zeros(Float64,Nrad)
  nash = zeros(Float64,Nrad)
  f0ash = zeros(Float64,Nrad,Nv)
  Tash = zeros(Float64,Nrad)

  for ir in 1:Nrad
    ir_set = ir
    logf0 = vec(log(abs(f0alpha[ir,:])))
    vt_temp = vt_temp_fac*sqrt(2.0*Ti[ir]*Tashfac/m_trace)
    Nv_use = indmin(abs(v-vt_temp))
#    coeffs = linear_fit(v[1:Nv_use].^2,logf0[1:Nv_use])
    coeffs = exp_fit(v[1:Nv_use].^2,abs(vec(f0alpha[ir,1:Nv_use])))
#    coeffs = exp_fit(v[1:Nv_use].^2,vec(f0alpha[ir,1:Nv_use]))
    Tash[ir] = abs(0.5*m_trace/coeffs[2])
    nash[ir] = coeffs[1]*(2.0*pi*Tash[ir]/m_trace)^1.5
    f0ash[ir,:] = coeffs[1]*exp(coeffs[2]*v.^2)
    ashfit[ir] = sqrt( sum( vec(f0alpha[ir,1:Nv_use]-f0ash[ir,1:Nv_use]).^2./(f0ash[ir,1:Nv_use].^2)))/Nv
  end

  if adj_postproc_temp 
    f0alpha, f0ash = rescale_ashtemp(f0alpha,nash,Tash,Ti)
  end
  Tash = Ti
 
  dfdr = Array(Float64,(Nrad,Nv))
  dfdv = Array(Float64,(Nrad,Nv))
  for iv in 1:Nv
    dfdr[1,iv] = (f0alpha[2,iv] - f0alpha[1,iv])/(rgrid[2]-rgrid[1])
    for ir in 2:Nrad-1
      dfdr[ir,iv] = (f0alpha[ir+1,iv] - f0alpha[ir-1,iv])/(rgrid[ir+1]-rgrid[ir-1])
    end
    dfdr[end,iv] = (f0alpha[end,iv] - f0alpha[end-1,iv])/(rgrid[end]-rgrid[end-1])
  end
  for ir in 1:Nrad
    dfdv[ir,1] = (f0alpha[ir,2] - f0alpha[ir,1])/(v[2]-v[1])
    for iv in 2:Nv-1
      dfdv[ir,iv] = (f0alpha[ir,iv+1] - f0alpha[ir,iv-1])/(v[iv+1]-v[iv-1])
    end
    dfdv[ir,end] = (f0alpha[ir,end] - f0alpha[ir,end-1])/(v[end]-v[end-1])
  end
  rflux = Array(Float64,(Nrad,Nv))
  vflux = Array(Float64,(Nrad,Nv))
  vflux_turb = Array(Float64,(Nrad,Nv))
  flux0 = Array(Float64,Nrad)
  hflux0 = Array(Float64,Nrad)
  pflux = Array(Float64,Nrad)
  hflux = Array(Float64,Nrad)
  Dr = Array(Float64,Nrad)
  chieff = Array(Float64,Nrad)
  for ir in 1:Nrad
    for iv in 1:Nv
      rflux[ir,iv] = -grad_rho[ir]*(Drr[ir,iv]*dfdr[ir,iv] + Drv[ir,iv]*dfdv[ir,iv])
      vflux[ir,iv] = -Dvr[ir,iv]*dfdr[ir,iv] - Dvv[ir,iv]*dfdv[ir,iv] - nu_s_v3[ir,iv]*f0alpha[ir,iv]/v[iv]^2 - 0.5*nu_par_v3[ir,iv]*dfdv[ir,iv]/v[iv]^2
      vflux_turb[ir,iv] = -Dvr[ir,iv]*dfdr[ir,iv] - Dvv[ir,iv]*dfdv[ir,iv] 
    end
    flux0[ir] = grad_rho[ir]*dot(d3v,-vec(Drv[ir,:].*dfdv[ir,:]))
    pflux[ir] = dot(d3v,vec(rflux[ir,:]))
    hflux[ir] = dot(d3v,0.5*m_trace*vec(rflux[ir,:]).*v.^2)
    Dr[ir] = (flux0[ir] - pflux[ir])/dot(d3v,vec(dfdr[ir,:]))
    hflux0[ir] = grad_rho[ir]*dot(0.5*m_trace*d3v.*v.^2,-vec(Drv[ir,:].*dfdv[ir,:]))
    chieff[ir] = (hflux0[ir] - hflux[ir])/dot(0.5*m_trace*d3v.*(v.^2),vec(dfdr[ir,:]))
  end

  nus_out = Array(Float64,(Nrad,Nv))
  nupar_out = Array(Float64,(Nrad,Nv))
  nupar_v2 = Array(Float64,(Nrad,Nv))
  Delta = Array(Float64,(Nrad,Nv))
  for ir in 1:Nrad
    nus_out[ir,:] = vec(nu_s_v3[ir,:])./(v.^3)
    nupar_out[ir,:] = vec(nu_par_v3[ir,:])./(v.^3)
    nupar_v2[ir,:] = 0.5*vec(nu_par_v3[ir,:])./v
    Delta = (Drv+Dvr).^2 - 4.0*Drr.*(nupar_v2+Dvv)
  end

  nhot = zeros(Float64,Nrad)
  nsd = zeros(Float64,Nrad)
  Tsd = zeros(Float64,Nrad)
  nalpha = zeros(Float64,Nrad)
  f0sd = zeros(Float64,Nrad,Nv)
  f0sdpure = zeros(Float64,Nrad,Nv)
  f0hot = zeros(Float64,Nrad,Nv)
  Talpha = zeros(Float64,Nrad)
  totheating = zeros(Float64,Nrad)
  collcreate = zeros(Float64,Nrad)
  ionheating = zeros(Float64,Nrad)
  elheating = zeros(Float64,Nrad)
  ashheating = zeros(Float64,Nrad)
  ihotheating = zeros(Float64,Nrad)
  ehotheating = zeros(Float64,Nrad)
  hotheating = zeros(Float64,Nrad)
  analheating = zeros(Float64,Nrad)
  sdiheating = zeros(Float64,Nrad)
  sdeheating = zeros(Float64,Nrad)
  totheatingv = zeros(Float64,(Nrad,Nv))
  iheatingv = zeros(Float64,(Nrad,Nv))
  eheatingv = zeros(Float64,(Nrad,Nv))
  sdiheatingv = zeros(Float64,(Nrad,Nv))
  sdeheatingv = zeros(Float64,(Nrad,Nv))
  sdpureiheatingv = zeros(Float64,(Nrad,Nv))
  sdpureeheatingv = zeros(Float64,(Nrad,Nv))
  energy = 0.5*m_trace*v.^2

  for ir in 1:Nrad
    nalpha[ir]=dot(vec(f0alpha[ir,:]),d3v)
#    f0sd[ir,:] = find_local_sd(ir,nalpha[ir])
#    f0sd[ir,:] = analytic_sd(ir,nalpha[ir],Ti[ir],false)
    f0sd[ir,:] = broad_sd(ir,nalpha[ir],Ti[ir],false)
#    f0sd[ir,:] = analytic_sd(ir,nalpha[ir],Tash[ir],false)
    f0sdpure[ir,:] = analytic_sd(ir,nash[ir],Tash[ir],true)
    nsd[ir]=dot(vec(f0sd[ir,:]),d3v)
    Tsd[ir]=dot(0.5*m_trace*v.^2.*vec(f0sd[ir,:]),d3v)
#    f0sd[ir,:] += f0ash[ir,:]
    # For ash-fit debugging purposes:

    f0hot[ir,:] = f0alpha[ir,:] - f0ash[ir,:]
    nhot[ir] = dot(d3v,vec(f0hot[ir,:]))
    energy = 0.5*m_trace*v.^2
    Talpha[ir] = dot(d3v,energy.*vec(f0alpha[ir,:]))/nalpha[ir]
    totheating[ir] = -dot(d3v.*energy,collop[:,:,ir]*vec(f0alpha[ir,:]))
    collcreate[ir] = -dot(d3v,collop[:,:,ir]*vec(f0alpha[ir,:]))
    ionheating[ir] = -dot(d3v.*energy,collop_ion[:,:,ir]*vec(f0alpha[ir,:]))
    elheating[ir] = -dot(d3v.*energy,collop_el[:,:,ir]*vec(f0alpha[ir,:]))
    ashheating[ir] = -dot(d3v.*energy,collop[:,:,ir]*vec(f0ash[ir,:]))
    hotheating[ir] = -dot(d3v.*energy,collop[:,:,ir]*vec(f0hot[ir,:]))
    ihotheating[ir] = -dot(d3v.*energy,collop_ion[:,:,ir]*vec(f0hot[ir,:]))
    ehotheating[ir] = -dot(d3v.*energy,collop_el[:,:,ir]*vec(f0hot[ir,:]))
#    analheating[ir],err = quadgk(heat_integrand,0.0,vmax)
    analheating[ir],err = hquadrature(heat_integrand,0.0,vmax)
    sdiheating[ir] = -dot(d3v.*energy,collop_ion[:,:,ir]*vec(f0sd[ir,:]))
    sdeheating[ir] = -dot(d3v.*energy,collop_el[:,:,ir]*vec(f0sd[ir,:]))

    totheatingv[ir,:] = -0.5*m_trace*4.0*pi*(v.^4).*(collop[:,:,ir]*vec(f0alpha[ir,:]))
    iheatingv[ir,:] = -0.5*m_trace*4.0*pi*(v.^4).*(collop_ion[:,:,ir]*vec(f0alpha[ir,:]))
    eheatingv[ir,:] = -0.5*m_trace*4.0*pi*(v.^4).*(collop_el[:,:,ir]*vec(f0alpha[ir,:]))
    sdiheatingv[ir,:] = -0.5*m_trace*4.0*pi*(v.^4).*(collop_ion[:,:,ir]*vec(f0sd[ir,:]))
    sdeheatingv[ir,:] = -0.5*m_trace*4.0*pi*(v.^4).*(collop_el[:,:,ir]*vec(f0sd[ir,:]))
    sdpureiheatingv[ir,:] = -0.5*m_trace*4.0*pi*(v.^4).*(collop_ion[:,:,ir]*vec(f0sdpure[ir,:]))
    sdpureeheatingv[ir,:] = -0.5*m_trace*4.0*pi*(v.^4).*(collop_el[:,:,ir]*vec(f0sdpure[ir,:]))
  end

  tot_integrated_heating = sum((rgrid[2]-rgrid[1])*Vprime.*totheating)
  ion_integrated_heating = sum((rgrid[2]-rgrid[1])*Vprime.*ionheating)
  #println("Integrated heating = ",integrated_heating)
  println("Integrated ion heating = ",ion_integrated_heating)

 
  dfashdr = Array(Float64,(Nrad,Nv))
  dfashdv = Array(Float64,(Nrad,Nv))
  for iv in 1:Nv
    dfashdr[1,iv] = (f0ash[2,iv] - f0ash[1,iv])/(rgrid[2]-rgrid[1])
    for ir in 2:Nrad-1
      dfashdr[ir,iv] = (f0ash[ir+1,iv] - f0ash[ir-1,iv])/(rgrid[ir+1]-rgrid[ir-1])
    end
    dfashdr[end,iv] = (f0ash[end,iv] - f0ash[end-1,iv])/(rgrid[end]-rgrid[end-1])
  end
  for ir in 1:Nrad
    dfashdv[ir,1] = (f0ash[ir,2] - f0ash[ir,1])/(v[2]-v[1])
    for iv in 2:Nv-1
      dfashdv[ir,iv] = (f0ash[ir,iv+1] - f0ash[ir,iv-1])/(v[iv+1]-v[iv-1])
    end
    dfashdv[ir,end] = (f0ash[ir,end] - f0ash[ir,end-1])/(v[end]-v[end-1])
  end

  rfluxash = Array(Float64,(Nrad,Nv))
  pfluxash = Array(Float64,Nrad)
  hfluxash = Array(Float64,Nrad)
  for ir in 1:Nrad
    for iv in 1:Nv
      rfluxash[ir,iv] = -grad_rho[ir]*(Drr[ir,iv]*dfashdr[ir,iv] + Drv[ir,iv]*dfashdv[ir,iv])
    end
    pfluxash[ir] = dot(d3v,vec(rfluxash[ir,:]))
    hfluxash[ir] = dot(d3v,0.5*m_trace*vec(rfluxash[ir,:]).*(v.^2))
  end

  heatfile = open("totheating.dat","a")
  writedlm(heatfile,[nedge,turbfac,tot_integrated_heating,ion_integrated_heating]')
  close(heatfile)

  sourcetot = Array(Float64,length(rgrid_in))
  for ir in 1:length(rgrid_in)
    sourcetot[ir] = dot(d3v,vec(source_in[ir,:]))
  end

  vref = sqrt(2.0*Tref/mref)
 
  # Save 2D plots
  plotnsave_2d(f0alpha,"f0alpha",L"$F_0$",true)
  plotnsave_2d(f0alpha,"f0lin",L"$F_0$",false)
  plotnsave_2d(f0sd,"f0sd",L"$F_s$",true)
  plotnsave_2d(f0ash,"f0ash",L"$F_\mathrm{ash}$",false)
  plotnsave_2d(f0hot,"f0hot",L"$F_\mathrm{hot}$",true)
  plotnsave_2d(Drr,"Drr",L"$D_{\psi \psi}$",false)
  plotnsave_2d(Drv,"Drv",L"$D_{\psi v}$",false)
  plotnsave_2d(Dvr,"Dvr",L"$D_{v \psi}$",false)
  plotnsave_2d(Dvv,"Dvv",L"$D_{v v}$",false)
  plotnsave_2d(rflux,"rflux",L"$\Gamma_\psi$",false)
  plotnsave_2d(rfluxash,"rfluxash","rfluxash",false)
  plotnsave_2d(vflux,"vflux",L"$\Gamma_v$",false)
  plotnsave_2d(source_local,"source",L"$S(r,v)$",false)
  plotnsave_2d(nus_out,"nus",L"$\nu_s$",false)
  plotnsave_2d(nupar_out,"nupar",L"$\nu_\|$",false)
  plotnsave_2d(nupar_v2,"nuparv2",L"$v^2 \nu_\| / 2 $",false)
  plotnsave_2d(Delta,"descriminant",L"$\Delta$",false)
  plotnsave_2d(totheatingv,"totheatingv","totheatingv",false)
  plotnsave_2d(iheatingv,"iheatingv","iheatingv",false)
  plotnsave_2d(eheatingv,"eheatingv","eheatingv",false)
  plotnsave_2d(sdiheatingv,"sdiheatingv","sdiheatingv",false)
  plotnsave_2d(sdeheatingv,"sdeheatingv","sdeheatingv",false)
  plotnsave_2d(sdpureiheatingv,"sdpureiheatingv","sdpureiheatingv",false)
  plotnsave_2d(sdpureeheatingv,"sdpureeheatingv","sdpureeheatingv",false)
  plotnsave_2d(dfdv,"dfdv","dfdv",false)
  plotnsave_2d(dfdr,"dfdr","dfdr",false)
 
  # Save functions of v
  plotnsave_v(v,"vgrid","vgrid",false)
  plotnsave_v(d3v,"d3v","d3v",false)
  plotnsave_v(sum(f0alpha,1),"losF0",L"Line of sight $f_0$",true)
  plotnsave_v(sum(f0sd,1),"losFsd",L"Line of sight $f_s$",true)

  # Save functions of r
  plotnsave_r(rgrid,"rgrid","rgrid",false)
  plotnsave_r(vcrit,"vcrit","vcrit",false)
  plotnsave_r(reaction_rate,"reactionrate",L"$\sigma$",false)
  plotnsave_r(nhot,"nhot",L"$n_\mathrm{hot}$",false)
  plotnsave_r(nash,"nash",L"$n_\mathrm{ash}$",false)
  plotnsave_r(nsd,"nsd",L"$n_\mathrm{sd}$",false)
  plotnsave_r(Tsd,"Tsd",L"$T_\mathrm{sd}$",false)
  plotnsave_r(nalpha,"nalpha",L"$n_\alpha$",false)
  plotnsave_r(Tash,"Tash",L"$T_\mathrm{ash}$",false)
  plotnsave_r(Talpha,"Talpha",L"$T_\alpha$",false)
  plotnsave_r(pflux,"pflux",L"$\Gamma_{p\alpha}$",false)
  plotnsave_r(flux0,"flux0",L"$\Gamma_0$",false)
  plotnsave_r(Dr,"Dr_eff",L"$D_\mathrm{eff}$",false)
  plotnsave_r(hflux0,"hflux0","hflux0",false)
  plotnsave_r(chieff,"chieff","chieff",false)
  plotnsave_r(hflux,"hflux",L"$q_{\alpha}$",false)
  plotnsave_r(hfluxash,"hfluxash",L"$q_{ash}$",false)
  plotnsave_r(pfluxash,"pfluxash",L"$\Gamma_{ash}$",false)
  plotnsave_r(totheating,"totheating","totheating",false)
  plotnsave_r(ionheating,"ionheating","totheating",false)
  plotnsave_r(elheating,"elheating","elheating",false)
  plotnsave_r(ashheating,"ashheating","ashheating",false)
  plotnsave_r(hotheating,"hotheating","hotheating",false)
  plotnsave_r(ihotheating,"ihotheating","ihotheating",false)
  plotnsave_r(ehotheating,"ehotheating","ehotheating",false)
  plotnsave_r(analheating,"analheating","analheating",false)
  plotnsave_r(ashfit,"ashfit","ashfit",false)
  plotnsave_r(taus,"taus",L"$\tau_s$",false)
  plotnsave_r(sdiheating,"sdiheating","sdiheating",false)
  plotnsave_r(sdeheating,"sdeheating","sdeheating",false)
  
  plotnsave_global(rgrid_in,"rgrid_global","rgrid",false)
  plotnsave_global(Te_in,"Te",L"$T_e$",false)
  plotnsave_global(Ti_in,"Ti",L"$T_i$",false)
  plotnsave_global(ne_in/1.0e20,"ne",L"$n_e$",false)
  plotnsave_global(sourcetot,"sourcetot",L"$\int S d^3\mathbf{v}$",false)
  plotnsave_global(surface_area_global,"area","Surface area",false)

  plotnsave_gs2(rgrid_gs2,"rgrid_gs2","rgrid_gs2",false)
  plotnsave_gs2(chii,"chii",L"$\chi_\mathrm{eff}$",false)
  plotnsave_gs2(phi2,"phi2",L"$|\phi^2|$",false)
  plotnsave_gs2(hflux_tot,"hflux_tot",L"$|\phi^2|$",false)
  plotnsave_gs2(rhostar,"rhostar",L"$\rho^*$",false)

  f0dat = Array(Float64,(Nrad*Nv,3))
  for ir in 1:Nrad
    for iv in 1:Nv
       f0dat[Nv*(ir-1)+iv,1] = rgrid[ir]/a
       f0dat[Nv*(ir-1)+iv,2] = v[iv] / valpha
       f0dat[Nv*(ir-1)+iv,3] = f0alpha[ir,iv]
    end
  end
     
  writedlm("f0.dat",f0dat)

  area_func = Spline1D(rgrid_in,surface_area_global)
  area_local = evaluate(area_func,rgrid[1])

  
  sourceenergy = Array{Float64}(Nrad)
  turbheating = Array{Float64}(Nrad)
  for ir in 1:Nrad
    sourceenergy[ir] = dot(vec(source_local[ir,:]),energy.*d3v)
    turbheating[ir] = m_trace*dot(vec(vflux_turb[ir,:]),v.*d3v)
  end
  cons_heatin = dot(fluxin.*energy,d3v)*Vprime[1]
  cons_source = dot(sourceenergy,Vprime)*(rgrid[2]-rgrid[1])
  cons_coll = dot(totheating,Vprime)*(rgrid[2]-rgrid[1])
  cons_coll_i = dot(ionheating,Vprime)*(rgrid[2]-rgrid[1])
  cons_coll_e = dot(elheating,Vprime)*(rgrid[2]-rgrid[1])
  cons_heatout = hflux[end]*Vprime[end]
  cons_turbheat = dot(turbheating,Vprime)*(rgrid[2]-rgrid[1])
  print("Heat flux in = ")
  println(cons_heatin)
  print("Energy generated by source = ")
  println(cons_source)
  print("Energy lost by collisions = ")
  println(cons_coll)
  print("Heating of ions = ")
  println(cons_coll_i)
  print("Heating of electrons = ")
  println(cons_coll_e)
  print("Heat flux out = ")
  println(cons_heatout)
  print("Turbulent heating of alphas = ")
  println(cons_turbheat)

  print("Conservation fraction = ")
  println(abs(( (cons_heatin + cons_source + cons_turbheat )-(cons_coll+cons_heatout ))/(cons_coll + cons_heatout)))
  
end

function rescale_ashtemp(f0alpha_in,nash,T_old,T_new)
   f0alpha_out = zeros(Float64,(Nrad,Nv))
   fash_out = zeros(Float64,(Nrad,Nv))
   for ir in 1:Nrad
      ashdist_old = exp(-m_trace*v.^2/(2.0*T_old[ir]))
      ashdist_old = ashdist_old*nash[ir]*(pi*2.0*T_old[ir]/m_trace)^(-1.5)
      ashdisk_old = ashdist_old'

      ashdist_new = exp(-m_trace*v.^2/(2.0*T_new[ir]))
      ashdist_new = ashdist_new*nash[ir]*(pi*2.0*T_new[ir]/m_trace)^(-1.5)
      ashdisk_new = ashdist_new'

      f0alpha_out[ir,:] = f0alpha_in[ir,:] - ashdist_old' + ashdist_new'
      fash_out[ir,:] = ashdist_new'
   end

   return f0alpha_out, fash_out
end

function plotnsave_v(data,varname,varlabel,logplot)
  global dir
  plotfile = dir*"/"*varname*".png"
  datfile = dir*"/"*varname*".dat"
  if logplot
    semilogy(v/valpha,vec(abs(data)))
  else
    plot(v/valpha,vec(data))
  end
  ylabel(varlabel)
  xlabel(L"$v / v_\alpha$")
  savefig(plotfile)
  cla()
  clf()
  
  writedlm(datfile,vec(data))
end

function plotnsave_global(data,varname,varlabel,logplot)
  global dir
  plotfile = dir*"/"*varname*".png"
  datfile = dir*"/"*varname*".dat"
  if logplot
    semilogy(rgrid_in/a,vec(abs(data)))
  else
    plot(rgrid_in/a,vec((data)))
  end
  ylabel(varlabel)
  xlabel(L"$r / a$")
  ylabel(varlabel)
  savefig(plotfile)
  cla()
  clf()
  
  writedlm(datfile,vec(data))
end

function plotnsave_gs2(data,varname,varlabel,logplot)
  global dir
  plotfile = dir*"/"*varname*".png"
  datfile = dir*"/"*varname*".dat"
  if logplot
    semilogy(rgrid_gs2/a,vec(abs(data)),"-o")
  else
    plot(rgrid_gs2/a,vec(data),"-o")
  end
  ylabel(varlabel)
  xlabel(L"$r / a$")
  ylabel(varlabel)
  savefig(plotfile)
  cla()
  clf()
  
  writedlm(datfile,vec(data))
end

function plotnsave_r(data,varname,varlabel,logplot)
  global dir
  datfile = dir*"/"*varname*".dat"
  if plot_output
  plotfile = dir*"/"*varname*".png"
  if logplot
    semilogy(rgrid/a,vec(abs(data)))
  else
    plot(rgrid/a,vec(data))
  end
  ylabel(varlabel)
  xlabel(L"$r / a$")
  ylabel(varlabel)
  savefig(plotfile)
  cla()
  clf()
  end
  
  writedlm(datfile,vec(data))
end

function save_density(f0_in)
   global f0save
   data = zeros(Nrad,2)
   data[:,1] = rgrid
   for ir in 1:Nrad
     data[ir,2] = dot(vec(f0in[ir,:]),d3v)/ne[ir]
   end
   writedlm("density.dat",data," ")
end
 
function heat_integrand(v)
  global ir_set, Tash, nash
    
    T_coll = [Ti[ir_set], Ti[ir_set], Te[ir_set]]
    n_coll = [DTmix*ne[ir_set], (1.0-DTmix)*ne[ir_set], ne[ir_set]]
    m_coll = [2.0*mp, 3.0*mp, me]
    Z_coll = [1.0,1.0,-1.0]
 
    nu_s_v3 = 0.0
    nu_par_v5 = 0.0
    for is in 1:3
      logLambda = lnLambda(m1=m_trace, m2=m_coll[is], Z1=Z_trace, Z2=Z_coll[is], T1=10.0*Ti[ir_set], T2=T_coll[is], n1=0.01*ne[ir_set], n2=n_coll[is])
      vts = sqrt(2.0*T_coll[is]/m_coll[is])
      
      nuhat_vta3 = n_coll[is]*Z_coll[is]^2*Z_trace^2*el^4*logLambda
      nuhat_vta3 = nuhat_vta3 / (4.0*pi*(ep0*m_trace)^2)

      nu_s_v3 = nu_s_v3 + nuhat_vta3*m_trace*G(v/vts)*v*v/T_coll[is]

      nu_par_v5 = nu_par_v5 + 2.0*nuhat_vta3*G(v/vts)*v^2
    end
    temp= - 4.0*pi*m_trace*v*(1.0 - 0.5m_trace*v^2/Tash[ir_set])*nash[ir_set]*(m_trace/(2.0*pi*Tash[ir_set]))^1.5*exp(-0.5*m_trace*v^2/Tash[ir_set])*(nu_s_v3 - 0.5*nu_par_v5*m_trace/Tash[ir_set])
    return temp
end
 
function plotnsave_2d(data,varname,varlabel,logplot)
  global dir

  datfile = dir*"/"*varname*".dat"
  if plot_output
  plotfile = dir*"/"*varname*".png"
  localfile = dir*"/"*varname*"_local.png"
  if logplot
    pcolormesh(rgrid/a,v/valpha,log(abs(data')))
  else
    pcolormesh(rgrid/a,v/valpha,data')
  end
  xlabel(L"$r / a$")
  ylabel(L"$v / v_\alpha$")
  title(varlabel)
  colorbar()
  savefig(plotfile)
  cla()
  clf()

  # Also save a plot local to particular radius
  if logplot
    semilogy(v/valpha,vec(abs(data[ir_sample,:])))
  else
    plot(v/valpha,vec(data[ir_sample,:]))
  end
  xlabel(L"$v / v_\alpha$")
  ylabel(varlabel)
  savefig(localfile)
  cla()
  clf()
  end
   
  writedlm(datfile,data)
end


end
