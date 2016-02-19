module GS2wrapper
using NetCDF
include("Species.jl"); using Species
include("input.jl")

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
      charge = ncread(filename,"z")[is] * qref
      dens = ncread(filename,"dens")[is] * nref
      temp = ncread(filename,"temp")[is] * Tref
      tprim = ncread(filename,"tprim")[is] * a / Tref
      fprim = ncread(filename,"fprim")[is] * a / nref
      push!(bulkspec, SpeciesData(mass=mass,q=charge,dens=dens,temp=temp,tprim=tprim,fprim=fprim,ispec=is) )
   end
 
   tracespec = SpeciesData[]
   for is in traceidx
      mass = ncread(filename,"mass")[is] * mref
      charge = ncread(filename,"z")[is] * qref
      dens = ncread(filename,"dens")[is] * nref
      temp = ncread(filename,"temp")[is] * Tref
      tprim = ncread(filename,"tprim")[is] * a / Tref
      fprim = ncread(filename,"fprim")[is] * a / nref
      push!(tracespec, SpeciesData(mass=mass,q=charge,dens=dens,temp=temp,tprim=tprim,fprim=fprim,ispec=is) )
      (tracespec[is].tprim == 0.0) || println("WARNING: trace species should have no temperature gradient to avoid singularity of solution")
      if is > 1
         (tracespec[is].mass == tracespec[is-1].mass) || error("Trace species masses must match")
         (tracespec[is].q == tracespec[is-1].q) || error("Trace species charges must match")
      end
   end
     
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
      egrid = vec(ncread(filename,"egrid")[:,spec.is]) * spec.temp
   end
   Nt = length(time)
   NE = length(egrid)

   rflux = zeros(Float64,NE)

   for ie in 1:NE
     rflux[ie] += try time_avg(vec(ncread(filename,"es_flux_e")[ie,spec.ispec,:]),time)
     rflux[ie] += try time_avg(vec(ncread(filename,"apar_flux_e")[ie,spec.ispec,:]),time)
     rflux[ie] += try time_avg(vec(ncread(filename,"bpar_flux_e")[ie,spec.ispec,:]),time)
   end

   return rflux * rhostar^2 * nref * vref
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
      egrid = vec(ncread(filename,"egrid")[:,spec.is]) * spec.temp
   end
   Nt = length(time)
   NE = length(egrid)
   vgrid = sqrt(2.0*egrid/spec.mass)

   vflux = zeros(Float64,NE)

   for ie in 1:NE
     vflux[ie] += try time_avg(vec(ncread(filename,"es_eflux")[ie,spec.ispec,:]),time)
     vflux[ie] += try time_avg(vec(ncread(filename,"apar_eflux")[ie,spec.ispec,:]),time)
     vflux[ie] += try time_avg(vec(ncread(filename,"bpar_eflux")[ie,spec.ispec,:]),time)
   end

   return vflux * rhostar^2 * nref * vref * Tref ./ (a * spec.mass * vgrid)
end

function get_egrid_from_GS2(filename::AbstractString,spec::SpeciesData)
   if agk
      egrid = vec(ncread(filename,"egrid")) * spec.temp
   else
      egrid = vec(ncread(filename,"egrid")[:,spec.is]) * spec.temp
   end

   return egrid
end

function get_df0dE_from_GS2(filename::AbstractString,spec::SpeciesData)
   if agk
      egrid = vec(ncread(filename,"egrid")) * spec.temp
      f0 = spec.dens * (spec.mass/(2.0*pi*spec.temp))^1.5 * exp(-egrid/spec.temp)
      df0dE = -f0/spec.temp
   else
      df0dE = vec(ncread(filename,"df0dE")[:,spec.is]) .*vec(ncread(filename,"f0")[:,spec.is]) * spec.dens
   end

   return df0dE
end


function time_avg(y::Array{Float64,1},t::Array{Float64,1})
   # Time to begin averaging
   tavg = t[1] + tavg_frac*(t[end]-t[1])
   
   # Find the closest index
   it = indmin(abs(t-tavg))

   return sum( (time[it+1:end]-time[it:end-1]).*y[it+1:end] ) / (time[end]-time[it])
end

end
