module Species
using Dierckx
using CurveFit
include("Input.jl")
include("Grids.jl"); using Grids
#include("Constants.jl")
#using .Constants

export SpeciesData

"A container for various plasma species data. All quantities must be specified."
type SpeciesData
   mass::Real
   q::Real
   dens::Real
   temp::Real
   fprim::Real
   tprim::Real
   ispec::Integer
end
SpeciesData(;mass=-1.0,q=-1.0,dens=-1.0,temp=-1.0,fprim=1.0,tprim=1.0,ispec=1) =SpeciesData(mass,q,dens,temp,fprim,tprim,ispec)

function calculate_diffcoeffs

end

function calculate_Drr

end

"""
Dvv = calculate_Dvv()

Calculates the transport coefficient Dvv based on energy fluxes for the trace species
"""
function calculate_Dvv(specs::Array{SpeciesData,1},vgrid_gs2::Array{Float64,1},vflux_gs2::Array{Float64,2},vgrid::Grid)
   Nspec, Nv_gs2 = size(vfluxes)

   # Interpolate fluxes to the t3core vgrid from the gs2 vgrid
   vflux_func = Spline1D(vgrid_gs2, flux_gs2, k=1,bc="extrapolate")
   vflux = evaluate(vflux_func,vgrid.pts)

   # Calculate Dvv and "H0" by finding linear fit of fluxes to df0/dv
   Dvv = zeros(Float64,Nv)
   H0 = zeros(Float64,Nv)
     
   return Dvv, H0
end



end

