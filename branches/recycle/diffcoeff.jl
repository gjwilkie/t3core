module diffcoeff
using input:Nrad,Nv, Nrad_gs2, tracespecs, tavg, rgrid_gs2, a, rhostar, mref, diffmodel, D0, brv,bvr, ir_sample, dilution_model, vflux_fac, m_trace, diff_D0, diff_v0, diff_power, change_diffmodel, constantD, turbfac, emrescale, Tashfac_in, spline_k, dilute_fac
using constants: valpha
using grids: v, rgrid, d3v
using matrix
using Dierckx
using species: filenames,temperature, density, mass, dens_gs2, temp_gs2, nref, Tref, Ti
using NetCDF
using PyPlot

#plt.style[:use]("myfig")

export Drv,Drr,Dvr,Dvv, init_diffCoeff, init_diffCoeff_zero, Drr_gs2, Drv_gs2, Dvr_gs2, Dvv_gs2, egrid_gs2, vgrid_gs2, Dnn, DnT, DTn, DTT, pflux0, hflux0, Dnn_gs2, DnT_gs2, DTn_gs2, DTT_gs2, pflux0_gs2, hflux0_gs2, init_diffCoeff_maxw, chii, phi2, hflux_tot

Drr = Float64[]
Drv = Float64[]
Dvr = Float64[]
Dvv = Float64[]
Drr_gs2 = Float64[]
Drv_gs2 = Float64[]
Dvr_gs2 = Float64[]
Dvv_gs2 = Float64[]
egrid_gs2 = Float64[]
vgrid_gs2 = Float64[]

Dnn_gs2 = Float64[]
DnT_gs2 = Float64[]
DTn_gs2 = Float64[]
DnT_gs2 = Float64[]
pflux0_gs2 = Float64[]
hflux0_gs2 = Float64[]
Dnn = Float64[]
DnT = Float64[]
DTn = Float64[]
DTT = Float64[]
pflux0 = Float64[]
hflux0 = Float64[]
chii = Float64[]
phi2 = Float64[]
hflux_tot = Float64[]

function init_diffCoeff_zero()
  global Drr, Drv, Dvr, Dvv
  Drr = zeros(Float64,(Nrad,Nv))
  Drv = zeros(Float64,(Nrad,Nv))
  Dvr = zeros(Float64,(Nrad,Nv))
  Dvv = zeros(Float64,(Nrad,Nv))
end

function init_diffCoeff()
  global Drr, Drv, Dvr, Dvv, Drr_gs2, Drv_gs2, Dvr_gs2, Dvv_gs2, egrid_gs2, vgrid_gs2, chii, phi2, hflux_tot

  ne_gs2 = length(ncread(filenames[1],"egrid")[:,1])

  Drr_gs2 = zeros(Nrad_gs2,Nv)
  Drv_gs2 = zeros(Nrad_gs2,Nv)
  Dvr_gs2 = zeros(Nrad_gs2,Nv)
  Dvv_gs2 = zeros(Nrad_gs2,Nv)
  Drr = zeros(Nrad,Nv)
  Drv = zeros(Nrad,Nv)
  Dvr = zeros(Nrad,Nv)
  Dvv = zeros(Nrad,Nv)

  egrid_gs2=Array(Float64,(ne_gs2))
  rflux_gs2=Array(Float64,(ne_gs2))
  rflux=Array(Float64,(Nv,length(tracespecs)))
  vflux_gs2=Array(Float64,(ne_gs2))
  vflux=Array(Float64,(Nv,length(tracespecs)))
  phi2 = Array(Float64,Nrad_gs2)
  hflux_tot = Array(Float64,Nrad_gs2)
  chii = Array(Float64,Nrad_gs2)
 
  ms = mass[tracespecs[2]]
  egrid = 0.5*ms*v.^2

  vref = sqrt(2.0*Tref/mref)

  apar_on = true
  bpar_on = true
  eflux_on = true

  for ir in 1:Nrad_gs2
    try(ncread(filenames[ir],"apar_flux_e"))
    catch
      apar_on = false
    end

    try(ncread(filenames[ir],"bpar_flux_e"))
    catch
      bpar_on = false
    end

    try(ncread(filenames[ir],"es_eflux"))
    catch
      eflux_on = false
    end
  end

  if !eflux_on && (diffmodel == 1)
    println("WARNING: Energy flux cannot be read from all nc files. Setting diffmodel to 2")
    change_diffmodel(2)
  end

  for ir in 1:Nrad_gs2
    fluxref = rhostar[ir]^2*nref[ir]*vref[ir]
    for is in 1:length(tracespecs)
      
      egrid_gs2[:]=ncread(filenames[ir],"egrid")[:,tracespecs[is]]*temp_gs2[ir,tracespecs[is]]
      vgrid_gs2=sqrt(2.0*egrid_gs2/mass[tracespecs[is]])

      time=ncread(filenames[ir],"t")[:]

      if (time[end]-time[1]) < tavg[ir]
        println("time[end] - time[1] = ",time[end]-time[1])
        println("tavg = ",tavg[ir])
        error("Must set tavg to be an appropriate window over which to time-average every GS2 simulation. This must be less than the total simulation time, which it isn't for radius $ir")
      end

      it = indmin(abs(time-(time[end]-tavg[ir])))

      # Time-averaged flux(E)
      for ie in 1:ne_gs2
        rflux_gs2[ie] =  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"es_flux_e")[ie,tracespecs[is],it:end-1]))/(time[end]-time[it])
        vflux_gs2[ie] =  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"es_eflux")[ie,tracespecs[is],it:end-1]))/(time[end]-time[it])
        if apar_on
          rflux_gs2[ie] +=  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"apar_flux_e")[ie,tracespecs[is],it:end-1]))/(time[end]-time[it])
          vflux_gs2[ie] +=  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"apar_eflux")[ie,tracespecs[is],it:end-1]))/(time[end]-time[it])
        end
        if bpar_on
          rflux_gs2[ie] +=  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"bpar_flux_e")[ie,tracespecs[is],it:end-1]))/(time[end]-time[it])
          vflux_gs2[ie] +=  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"bpar_eflux")[ie,tracespecs[is],it:end-1]))/(time[end]-time[it])
        end
        vflux_gs2[ie] = vflux_gs2[ie]*Tref[ir]/(a*ms*vgrid_gs2[ie])
      end

      phi2[ir] =  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"phi2")[it:end-1]))/(time[end]-time[it])
      hflux_tot[ir] =  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"hflux_tot")[it:end-1]))/(time[end]-time[it])
      hfluxi =  sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"es_heat_flux")[1,it:end-1]))/(time[end]-time[it])
      if apar_on 
        hfluxi += sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"apar_heat_flux")[1,it:end-1]))/(time[end]-time[it])
      end
      if bpar_on 
        hfluxi += sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"bpar_heat_flux")[1,it:end-1]))/(time[end]-time[it])
      end

      chii[ir] = hfluxi*rhostar[ir]^2*a*vref[ir]/(ncread(filenames[ir],"tprim")[1]*ncread(filenames[ir],"dens")[1]*ncread(filenames[ir],"temp")[1])

      rflux_func = Spline1D(vgrid_gs2,rflux_gs2,k=1,bc="nearest")
      vflux_func = Spline1D(vgrid_gs2,vflux_gs2,k=1,bc="nearest")

      rflux[:,is]=evaluate(rflux_func,v)
      vflux[:,is]=evaluate(vflux_func,v)

    end

    # Generate diffusion coefficients on GS2-r-grid
    temp1 = temp_gs2[ir,tracespecs[1]]
    temp2 = temp_gs2[ir,tracespecs[2]]
    dens1 = dens_gs2[ir,tracespecs[1]]
    dens2 = dens_gs2[ir,tracespecs[2]]
    tprim1 = ncread(filenames[ir],"tprim")[tracespecs[1]]/a
    tprim2 = ncread(filenames[ir],"tprim")[tracespecs[2]]/a
    fprim1 = ncread(filenames[ir],"fprim")[tracespecs[1]]/a
    fprim2 = ncread(filenames[ir],"fprim")[tracespecs[2]]/a

    if (abs(fprim1 - fprim2) < 10.0*eps(Float64)) && (abs(tprim1 - tprim2) < 10.0*eps(Float64)) && (abs(temp1 - temp2) < 10.0*eps(Float64)) 
      error("The two distributions used for diffusion model 1 must be different")
    end 
  
    gradf1 = fprim1 + ( 0.5*(ms*v.^2/temp1) - 1.5)*tprim1
    gradf2 = fprim2 + ( 0.5*(ms*v.^2/temp2) - 1.5)*tprim2

    d = ms*v.*((gradf1/temp2) - (gradf2/temp1))

    Drr_gs2[ir,:] = (ms/temp2)*v.*(rflux[:,1]/dens1) - (ms/temp1)*v.*(rflux[:,2]/dens2)
    Drv_gs2[ir,:] = -gradf2.*(rflux[:,1]/dens1) + gradf1.*(rflux[:,2]/dens2)
    Drr_gs2[ir,:] = fluxref*Drr_gs2[ir,:]./(d')
    Drv_gs2[ir,:] = fluxref*Drv_gs2[ir,:]./(d')

    Dvr_gs2[ir,:] = (ms/temp2)*v.*(vflux[:,1]/dens1) - (ms/temp1)*v.*(vflux[:,2]/dens2)
    Dvv_gs2[ir,:] = -gradf2.*(vflux[:,1]/dens1) + gradf1.*(vflux[:,2]/dens2)
    Dvr_gs2[ir,:] = fluxref*Dvr_gs2[ir,:]./(d')
    Dvv_gs2[ir,:] = fluxref*Dvv_gs2[ir,:]./(d')

  end

  for iv in 1:Nv
    Drr_func = Spline1D(rgrid_gs2, Drr_gs2[:,iv],k=spline_k)
    Drr[:,iv] = evaluate(Drr_func,rgrid)
    if (diffmodel < 3)
      Drv_func = Spline1D(rgrid_gs2, Drv_gs2[:,iv],k=spline_k)
      Drv[:,iv] = evaluate(Drv_func,rgrid)
      if (diffmodel == 1)
        Dvv_func = Spline1D(rgrid_gs2, Dvv_gs2[:,iv],k=spline_k)
        Dvr_func = Spline1D(rgrid_gs2, Dvr_gs2[:,iv],k=spline_k)
        Dvv[:,iv] = evaluate(Dvv_func,rgrid)
        Dvr[:,iv] = evaluate(Dvr_func,rgrid)
      end
    end
  end

  if  ((dilution_model == 1 ) || (dilution_model ==3)) & isfile("density.dat")
    data = readdlm("density.dat",' ')
    rgrid_d = data[:,1]
    dens_d = data[:,2]
    dens_func = Spline1D(rgrid_d,dens_d,k=spline_k)
    nhat = evaluate(dens_func,rgrid)
    for iv in 1:Nv
      Drr[:,iv] = Drr[:,iv] .* (1.0 - dilute_fac*nhat)
      Drv[:,iv] = Drv[:,iv] .* (1.0 - dilute_fac*nhat)
      Dvr[:,iv] = Dvr[:,iv] .* (1.0 - dilute_fac*nhat)
      Dvv[:,iv] = Dvv[:,iv] .* (1.0 - dilute_fac*nhat)
    end
  end

  if diffmodel == -1
    for ir in 1:Nrad
      vtash = sqrt(2.0*Ti[ir]*Tashfac_in/m_trace)
      rhostar_func = Spline1D(rgrid_gs2,rhostar.*sqrt(phi2),k=spline_k)
      rhostar_local = evaluate(rhostar_func,rgrid[ir])
      for iv in 1:Nv
        if v[iv] < diff_v0*vtash
          Drr[ir,iv] = diff_D0*rhostar_local^2*vtash*a
        else
          Drr[ir,iv] = diff_D0*rhostar_local^2*vtash*a*(v[iv]/(diff_v0*vtash))^diff_power
        end
      end
    end
    Drv[:,:] = 0.0
    Dvr[:,:] = 0.0
    Dvv[:,:] = 0.0
  elseif diffmodel == -2
    Drr[:,:] = constantD
    Drv[:,:] = 0.0
    Dvr[:,:] = 0.0
    Dvv[:,:] = 0.0
  end

  Drr = Drr*turbfac
  Drv = Drv*turbfac
  Dvr = Dvr*turbfac
  Dvv = Dvv*turbfac

  if emrescale
    for ir in 1:Nrad
      vtash = sqrt(2.0*Ti[ir]*Tashfac_in/m_trace)
      for iv in 1:Nv
        if v[iv] > diff_v0*vtash
          Drr[ir,iv] = Drr[ir,iv]*(v[iv]^2/(diff_v0*vtash)^2)
          Drv[ir,iv] = Drv[ir,iv]*(v[iv]^2/(diff_v0*vtash)^2)
          Dvr[ir,iv] = Dvr[ir,iv]*(v[iv]^2/(diff_v0*vtash)^2)
          Dvv[ir,iv] = Dvv[ir,iv]*(v[iv]^2/(diff_v0*vtash)^2)
        end
      end
    end
  end

end

function init_diffCoeff_maxw()
  global Dnn, DnT, DTn, DTT, pflux0, hflux0, Dnn_gs2, DnT_gs2, DTn_gs2, DTT_gs2, pflux0_gs2, hflux0_gs2

  Dnn_gs2 = zeros(Nrad_gs2)
  DnT_gs2 = zeros(Nrad_gs2)
  DTn_gs2 = zeros(Nrad_gs2)
  DTT_gs2 = zeros(Nrad_gs2)
  pflux0_gs2 = zeros(Nrad_gs2)
  hflux0_gs2 = zeros(Nrad_gs2)
  Dnn = zeros(Nrad)
  DnT = zeros(Nrad)
  DTn = zeros(Nrad)
  DTT = zeros(Nrad)
  pflux0 = zeros(Nrad)
  hflux0 = zeros(Nrad)

  pflux=Array(Float64,(length(tracespecs)))
  hflux=Array(Float64,(length(tracespecs)))
 
  ms = mass[tracespecs[2]]

  vref = sqrt(2.0*Tref/mref)

  for ir in 1:Nrad_gs2
    fluxref = rhostar[ir]^2*nref[ir]*vref[ir]
    for is in 1:length(tracespecs)
      
      time=ncread(filenames[ir],"t")[:]

      if (time[end]-time[1]) < tavg[ir]
        println("time[end] - time[1] = ",time[end]-time[1])
        println("tavg = ",tavg[ir])
        error("Must set tavg to be an appropriate window over which to time-average every GS2 simulation. This must be less than the total simulation time, which it isn't for radius $ir")
      end

      it = indmin(abs(time-(time[end]-tavg[ir])))

      # Time-averaged fluxes
      pflux[is] = fluxref* sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"es_part_flux")[tracespecs[is],it:end-1]))/(time[end]-time[it])
      hflux[is] = fluxref*Tref[ir]* sum( (time[it+1:end]-time[it:end-1]).*vec(ncread(filenames[ir],"es_heat_flux")[tracespecs[is],it:end-1]))/(time[end]-time[it])
      
    end

    # Generate diffusion coefficients on GS2-r-grid
    Dmatrix = Array(Float64,(3,3))
    for is in 1:3
       Dmatrix[is,1] = (nref[ir]/a)*ncread(filenames[ir],"dens")[tracespecs[is]] * ncread(filenames[ir],"fprim")[tracespecs[is]]
       Dmatrix[is,2] = (Tref[ir]/a)*ncread(filenames[ir],"temp")[tracespecs[is]] * ncread(filenames[ir],"tprim")[tracespecs[is]]
       Dmatrix[is,3] = 1.0
    end

    pvec = Dmatrix\pflux
    hvec = Dmatrix\hflux
    Dnn_gs2[ir] = pvec[1]
    DnT_gs2[ir] = pvec[2]
    pflux0_gs2[ir] = pvec[3]
    DTn_gs2[ir] = hvec[1]
    DTT_gs2[ir] = hvec[2]
    hflux0_gs2[ir] = hvec[3]
    
  end

  Dnn_func = Spline1D(rgrid_gs2, Dnn_gs2)
  Dnn[:] = evaluate(Dnn_func,rgrid)

  pflux0_func = Spline1D(rgrid_gs2, pflux0_gs2)
  pflux0[:] = evaluate(pflux0_func,rgrid)

  hflux0_func = Spline1D(rgrid_gs2, hflux0_gs2)
  hflux0[:] = evaluate(hflux0_func,rgrid)

  if (diffmodel < 3)
    DnT_func = Spline1D(rgrid_gs2, DnT_gs2[:])
    DnT[:] = evaluate(DnT_func,rgrid)
    if (diffmodel == 1)
      DTT_func = Spline1D(rgrid_gs2, DTT_gs2[:])
      DTn_func = Spline1D(rgrid_gs2, DTn_gs2[:])
      DTT[:] = evaluate(DTT_func,rgrid)
      DTn[:] = evaluate(DTn_func,rgrid)*vflux_fac
    end
  end

  if ( (dilution_model == 1 ) || (dilution_model == 3)) & isfile("density.dat")
    data = readdlm("density.dat",' ')
    rgrid_d = data[:,1]
    dens_d = data[:,2]
    dens_func = Spline1D(rgrid_d,dens_d)
    nhat = evaluate(dens_func,rgrid)
    for iv in 1:Nv
      Dnn[:] = Dnn[:] .* (1.0 - nhat)
      DnT[:] = DnT[:] .* (1.0 - nhat)
      DTn[:] = DTn[:] .* (1.0 - nhat)
      DTT[:] = DTT[:] .* (1.0 - nhat)
      pflux0[:] = pflux0[:] .* (1.0 - nhat)
      hflux0[:] = hflux0[:] .* (1.0 - nhat)
    end
  end


end

function init_diffcoeff_for_analytic_test()
  global Drv, Drr, Dvv, Dvr

  Drr = zeros(Nrad,Nv)
  Drv = zeros(Nrad,Nv)
  Dvr = zeros(Nrad,Nv)
  Dvv = zeros(Nrad,Nv)

  vhat = v/valpha
  rhat = rgrid/a

  for iv in 1:Nv
    for ir in 1:Nrad
      Drr[ir,iv] = D0*(3.0-rhat[ir])
      Drv[ir,iv] = D0*brv*exp(-2.0*vhat[iv])*valpha/a
      Dvr[ir,iv] = D0*bvr*exp(-2.0*vhat[iv])*valpha/a
      Dvv[ir,iv] = D0*0.25*exp(-vhat[iv])*valpha^2/a^2
    end
  end
  
end
end

