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

function calculate_diffcoeffs

end

function calculate_Drr

end

"""
Dvv = calculate_Dvv()

Calculates the transport coefficient Dvv based on energy fluxes for the trace species
"""
function calculate_Dvv(specs::Array{SpeciesData,1})

end

"""
rflux = get_rflux_from_GS2(filename::AbstractString,specs::SpeciesData)

Reads from the gs2 NetCDF output file filename and uses the defined trace species to return time-averaged radial flux
"""
function get_rflux_from_GS2(filename::AbstractString,spec::SpeciesData)
   time = vec(ncread(filename,"t"))
   egrid = vec(ncread(filename,"egrid"))
   Nt = length(time)
   NE = length(egrid)

   rflux = zeros(Float64,NE)

   for ie in 1:NE
     rflux[ie] += try time_avg(vec(ncread(filename,"es_flux_e")[ie,spec.ispec,:]),time)
     rflux[ie] += try time_avg(vec(ncread(filename,"apar_flux_e")[ie,spec.ispec,:]),time)
     rflux[ie] += try time_avg(vec(ncread(filename,"bpar_flux_e")[ie,spec.ispec,:]),time)
   end

   return rflux
end


"""
vflux = get_rflux_from_GS2(filename::AbstractString,specs::SpeciesData)

Reads from the gs2 NetCDF output file filename and uses the defined trace species to return time-averaged energy flux
"""
function get_rflux_from_GS2(filename::AbstractString,spec::SpeciesData)
   time = vec(ncread(filename,"t"))
   egrid = vec(ncread(filename,"egrid"))
   Nt = length(time)
   NE = length(egrid)

   vflux = zeros(Float64,NE)

   for ie in 1:NE
     vflux[ie] += try time_avg(vec(ncread(filename,"es_eflux")[ie,spec.ispec,:]),time)
     vflux[ie] += try time_avg(vec(ncread(filename,"apar_eflux")[ie,spec.ispec,:]),time)
     vflux[ie] += try time_avg(vec(ncread(filename,"bpar_eflux")[ie,spec.ispec,:]),time)
   end

   return vflux
end


function time_avg(y::Array{Float64,1},t::Array{Float64,1})
   # Time to begin averaging
   tavg = t[1] + tavg_frac*(t[end]-t[1])
   
   # Find the closest index
   it = indmin(abs(t-tavg))

   return sum( (time[it+1:end]-time[it:end-1]).*y[it+1:end] ) / (time[end]-time[it])
end

end
