module Species
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

end

