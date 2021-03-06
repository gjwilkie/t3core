module Species
include("Input.jl")
using Dierckx
using NetCDF
using CurveFit
using Grids

export SpeciesData, calculate_Dvv, Dvv_model, get_test_species, get_species_from_GS2, get_vgrid_from_GS2

"A container for various plasma species data. All quantities must be specified."
type SpeciesData
   mass::Float64
   q::Float64
   dens::Float64
   temp::Float64
   fprim::Float64
   tprim::Float64
   ispec::Int64
end
SpeciesData(;mass=-1.0,q=-1.0,dens=-1.0,temp=-1.0,fprim=1.0,tprim=1.0,ispec=1) =SpeciesData(mass,q,dens,temp,fprim,tprim,ispec)

function calculate_diffcoeffs

end

function calculate_Drr

end

"""
Calculates the transport coefficient Dvv based on energy fluxes for the trace species
"""
function calculate_Dvv(filename::AbstractString,tracespecs::Array{SpeciesData,1},vgrid::Array{Float64,1})
   Nspec = length(tracespecs)

   rflux = zeros(Float64,(Nspec,Nv))
   vflux = zeros(Float64,(Nspec,Nv))
   df0dr = zeros(Float64,(Nspec,Nv))
   df0dv = zeros(Float64,(Nspec,Nv))

   for is in 1:Nspec
      vgrid_gs2 = get_vgrid_from_GS2(filename,tracespecs[is])
      rflux_gs2 = get_rflux_from_GS2(filename,tracespecs[is]) 
      vflux_gs2 = get_vflux_from_GS2(filename,tracespecs[is]) 

      # Interpolate fluxes to the t3core vgrid from the gs2 vgrid
      rflux_func = Spline1D(vgrid_gs2, rflux_gs2, k=1,bc="extrapolate")
      rflux[is,:] = evaluate(rflux_func,vgrid)
      vflux_func = Spline1D(vgrid_gs2, vflux_gs2, k=1,bc="extrapolate")
      vflux[is,:] = evaluate(vflux_func,vgrid)

      df0dv_gs2 =  get_df0dv_from_GS2(filename,tracespecs[is])
      df0dv_func = Spline1D(vgrid_gs2, df0dv_gs2, k=1,bc="extrapolate")
      df0dv[is,:] = evaluate(df0dv_func,vgrid)

      df0dr_gs2 =  get_df0dr_from_GS2(filename,tracespecs[is])
      df0dr_func = Spline1D(vgrid_gs2, df0dr_gs2, k=1,bc="extrapolate")
      df0dr[is,:] = evaluate(df0dr_func,vgrid)

   end

   # Calculate Dvv and "H0" by finding linear fit of fluxes to df0/dv
   H0 = zeros(Float64,Nv)
   Dvv = zeros(Float64,Nv)

   if false

   fprim1 = tracespecs[1].fprim
   fprim2 = tracespecs[2].fprim
   tprim1 = tracespecs[1].tprim
   tprim2 = tracespecs[2].tprim
   temp1 = tracespecs[1].temp
   temp2 = tracespecs[2].temp
   dens1 = tracespecs[1].dens
   dens2 = tracespecs[2].dens
   ms = tracespecs[1].mass
   Drr = zeros(Nv)
   Drv = zeros(Nv)
   Dvr = zeros(Nv)
   Dvv = zeros(Nv)
   for iv in 1:Nv
      v = vgrid[iv]

      gradf1 = fprim1 + ( 0.5*(ms*v.^2/temp1) - 1.5)*tprim1
      gradf2 = fprim2 + ( 0.5*(ms*v.^2/temp2) - 1.5)*tprim2

      d = ms*v*((gradf1/temp2) - (gradf2/temp1))

      Drr[iv] = (ms/temp2)*v*(rflux[1,iv]/dens1) - (ms/temp1)*v.*(rflux[2,iv]/dens2)
      Drr[iv] = Drr[iv]/d
      Drv[iv] = -gradf2.*(rflux[1,iv]/dens1) + gradf1.*(rflux[2,iv]/dens2)
      Drv[iv] = Drv[iv]/d

      Dvr[iv] = (ms/temp2)*v*(vflux[1,iv]/dens1) - (ms/temp1)*v.*(vflux[2,iv]/dens2)
      Dvr[iv] = Dvr[iv]/d
      Dvv[iv] = -gradf2.*(vflux[1,iv]/dens1) + gradf1.*(vflux[2,iv]/dens2)
      Dvv[iv] = Dvv[iv]/d

   end
end

   dummy1, dummy2, dummy3, Dvv = getDiffCoeffs(df0dr,df0dv,rflux,vflux)
#   for iv in 1:Nv
#      a,b = linear_fit(df0dv[:,iv],vflux[:,iv])
#      H0[iv] = a
#      Dvv[iv] = -b
#   end
     
   return turb_rescale*H0, turb_rescale*Dvv
end

function getDiffCoeffs(df0dr::Array{Float64,2},df0dv::Array{Float64,2},rflux::Array{Float64,2},vflux::Array{Float64,2})
   Nv = length(df0dr[1,:])
   Drr = zeros(Float64,Nv)
   Drv = zeros(Float64,Nv)
   Dvr = zeros(Float64,Nv)
   Dvv = zeros(Float64,Nv)

   for iv in 1:Nv

      Drr[iv] = (df0dv[1,iv]*rflux[2,iv] - df0dv[2,iv]*rflux[1,iv])/(df0dv[2,iv]*df0dr[1,iv]-df0dv[1,iv]*df0dr[2,iv])
      Drv[iv] = -(rflux[1,iv]+Drr[iv]*df0dr[1,iv])/df0dv[1,iv]
      Dvr[iv] = (df0dv[1,iv]*vflux[2,iv] - df0dv[2,iv]*vflux[1,iv])/(df0dv[2,iv]*df0dr[1,iv]-df0dv[1,iv]*df0dr[2,iv])
      Dvv[iv] = -(vflux[1,iv]+Dvr[iv]*df0dr[1,iv])/df0dv[1,iv]
   end

#   dv = vmax / Nv
#   vgrid = collect( linspace(0.5*dv,vmax - 0.5*dv,Nv) )
#   semilogy(vgrid,abs(Drv))
#   semilogy(vgrid,abs(Dvr))
#   savefig("DrvVSDvr.png")
#   error("Stopping")
   println(Drr)
   return Drr, Drv, Dvr, Dvv
end

function Dvv_model(vgrid)

   v0 = vref*2
   D0 = 1.0 * (2.0*Tref/mref)/a^2
#   D0 = 2.0e12
#   v0 = 2.0e6

   H0 = zeros(Float64,Nv)
   Dvv = zeros(Float64,Nv)

   idx = indmin(abs(vgrid-v0))
   
   Dvv[1:idx-1] = D0

   # Fit to electrostatic fusion alphas:
   #Dvv[idx:end] = D0*(vgrid[idx:end]/v0).^(-4)
   # Extrapolation to active dBpar
   Dvv[idx:end] = D0*(vgrid[idx:end]/v0).^(-4)

#   H0 = Dvv .* (vgrid).^(-1)

   return turb_rescale*H0, turb_rescale*Dvv
   
end

function get_test_species()
   
   bulkspec = SpeciesData[]
      # ions
      mass = mref
      charge = qref
      dens = nref
      temp = Tref
      tprim = a / Tref
      fprim = a / nref
      is = 1
      push!(bulkspec, SpeciesData(mass=mass,q=charge,dens=dens,temp=temp,tprim=tprim,fprim=fprim,ispec=is) )

      # electrons
      mass = me
      charge = -qref
      dens = nref
      temp = Tref
      tprim = a / Tref
      fprim = a / nref
      is = 2
      push!(bulkspec, SpeciesData(mass=mass,q=charge,dens=dens,temp=temp,tprim=tprim,fprim=fprim,ispec=is) )
 
   tracespec = SpeciesData[]
   push!(tracespec, bulkspec[1])

   return bulkspec, tracespec
end

"""
bulkspecs, tracespecs = get_species_from_GS2(filename::AbstractString)

Reads the netcdf file filename and returns arrays containing the species data for bulk species and trace species (using the indices from Input.jl"
"""
function get_species_from_GS2(filename::AbstractString)
   nbulk = length(bulkidx) 
   ntrace = length(traceidx) 

   bulkspec = SpeciesData[]
   for is in bulkidx
      mass = ncread(filename,"mass")[is] * mref
      charge = ncread(filename,"charge")[is] * qref
      dens = ncread(filename,"dens")[is] * nref
      temp = ncread(filename,"temp")[is] * Tref
      tprim = ncread(filename,"tprim")[is]  / a
      fprim = ncread(filename,"fprim")[is]  / a
      push!(bulkspec, SpeciesData(mass=mass,q=charge,dens=dens,temp=temp,tprim=tprim,fprim=fprim,ispec=is) )
   end
 
   tracespec = SpeciesData[]
   for is in traceidx
      mass = ncread(filename,"mass")[is] * mref
      charge = ncread(filename,"charge")[is] * qref
      dens = ncread(filename,"dens")[is] * nref
      temp = ncread(filename,"temp")[is] * Tref
      tprim = ncread(filename,"tprim")[is] / a
      fprim = ncread(filename,"fprim")[is] / a
      push!(tracespec, SpeciesData(mass=mass,q=charge,dens=dens,temp=temp,tprim=tprim,fprim=fprim,ispec=is) )
      (tprim == 0.0) || println("WARNING: trace species should have no temperature gradient to avoid singularity of solution")
#      if is > 1
#         (mass == tracespec[is-1].mass) || error("Trace species masses must match")
#         (tracespec[is].q == tracespec[is-1].q) || error("Trace species charges must match")
#      end
   end

#   countnz(tracespec.mass-tracespec[1].mass) == 0 || error("Trace species masses must match")
#   countnz(tracespec.q-tracespec[1].q) == 0 || error("Trace species masses must match")
   length(tracespec) == 2 || error("Must use exactly two tracespecies for calculation to work")

#   print("bulkspec = ")
#   println(bulkspec)
#   print("tracespec = ")
#   println(tracespec)
     
   return bulkspec, tracespec
end

"""
rflux = get_rflux_from_GS2(filename::AbstractString,specs::SpeciesData)

Reads from the gs2 NetCDF output file filename and uses the defined trace species to return time-averaged radial flux
"""
function get_rflux_from_GS2(filename::AbstractString,spec::SpeciesData)
   time = vec(ncread(filename,"t"))
   if agk
      egrid = vec(ncread(filename,"egrid")) * spec.temp
   else
      egrid = vec(ncread(filename,"egrid")[:,spec.ispec]) * spec.temp
   end
   Nt = length(time)
   NE = length(egrid)

   rflux = zeros(Float64,NE)

   for ie in 1:NE
      if agk
         try rflux[ie] += time_avg(vec(ncread(filename,"es_rflux")[ie,spec.ispec,:]),time) end 
         try rflux[ie] += time_avg(vec(ncread(filename,"apar_rflux")[ie,spec.ispec,:]),time) end
         try rflux[ie] += time_avg(vec(ncread(filename,"bpar_rflux")[ie,spec.ispec,:]),time) end
      else
         try rflux[ie] += time_avg(vec(ncread(filename,"es_flux_e")[ie,spec.ispec,:]),time) end 
         try rflux[ie] += time_avg(vec(ncread(filename,"apar_flux_e")[ie,spec.ispec,:]),time) end
         try rflux[ie] += time_avg(vec(ncread(filename,"bpar_flux_e")[ie,spec.ispec,:]),time) end
      end
   end

   return rflux * rhostar^2 * vref * nref/spec.dens
end



"""
vflux = get_rflux_from_GS2(filename::AbstractString,specs::SpeciesData)

Reads from the gs2 NetCDF output file filename and uses the defined trace species to return time-averaged energy flux
"""
function get_vflux_from_GS2(filename::AbstractString,spec::SpeciesData)
   time = vec(ncread(filename,"t"))
   if agk
      egrid = vec(ncread(filename,"egrid")) * spec.temp
   else
      egrid = vec(ncread(filename,"egrid")[:,spec.ispec]) * spec.temp
   end
   Nt = length(time)
   NE = length(egrid)
   vgrid = sqrt(2.0*egrid/spec.mass)

   vflux = zeros(Float64,NE)

   for ie in 1:NE
      if agk
         try vflux[ie] += time_avg(vec(ncread(filename,"es_eflux")[ie,spec.ispec,:]),time) end 
         try vflux[ie] += time_avg(vec(ncread(filename,"apar_eflux")[ie,spec.ispec,:]),time) end
         try vflux[ie] += time_avg(vec(ncread(filename,"bpar_eflux")[ie,spec.ispec,:]),time) end
      else
         try vflux[ie] += time_avg(vec(ncread(filename,"es_eflux")[ie,spec.ispec,:]),time) end 
         try vflux[ie] += time_avg(vec(ncread(filename,"apar_eflux")[ie,spec.ispec,:]),time) end
         try vflux[ie] += time_avg(vec(ncread(filename,"bpar_eflux")[ie,spec.ispec,:]),time) end
      end
   end

   return vflux * rhostar^2 * vref * nref * Tref ./ (a * spec.mass * vgrid * spec.dens)
end

function get_vgrid_from_GS2(filename::AbstractString,spec::SpeciesData)
   if agk
      egrid = vec(ncread(filename,"egrid")) * spec.temp
   else
      egrid = vec(ncread(filename,"egrid")[:,spec.ispec]) * spec.temp
   end

   vgrid = sqrt(2.0*egrid / spec.mass)

   return vgrid
end

function get_df0dr_from_GS2(filename::AbstractString,spec::SpeciesData)
   vts = sqrt(2.0*spec.temp/spec.mass)
   if agk
      egrid = vec(ncread(filename,"egrid")) * spec.temp
      df0dr = -(ncread(filename,"fprim")[spec.ispec] + ncread(filename,"tprim")[spec.ispec] *(egrid/spec.temp - 1.5)) / a
   else
      df0dr = vec(ncread(filename,"f0prim")[:,spec.ispec]) / a
   end

   return df0dr
end


function get_df0dv_from_GS2(filename::AbstractString,spec::SpeciesData)
   if agk
      df0dE = -1.0/spec.temp
   else
      df0dE = vec(ncread(filename,"df0dE")[:,spec.ispec]) / Tref 
   end

   vgrid = get_vgrid_from_GS2(filename,spec)
   return df0dE.*(spec.mass*vgrid)
end


function time_avg(y::Array{Float64,1},t::Array{Float64,1})
   # Time to begin averaging
   tavg = t[1] + tavg_frac*(t[end]-t[1])
   
   # Find the closest index
   it = indmin(abs(t-tavg))

   return sum( (t[it+1:end]-t[it:end-1]).*y[it+1:end] ) / (t[end]-t[it])
end

end

