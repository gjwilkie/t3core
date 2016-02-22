module Species
using Dierckx
using CurveFit
include("Grids.jl"); using .Grids
include("Constants.jl")
include("Input.jl")

export SpeciesData, calculate_Dvv, Dvv_model, get_test_species

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
Calculates the transport coefficient Dvv based on energy fluxes for the trace species
"""
function calculate_Dvv(filename::AbstractString,tracespecs::Array{SpeciesData,1},vgrid::Array{Float64,1})
   Nspec = length(tracespecs)

   vflux = zeros(Float64,(Nspec,Nv))
   df0dv = zeros(Float64,(Nspec,Nv))

   for is in 1:Nspec
      vgrid_gs2 = get_vgrid_from_GS2(filename,tracespecs[is])
      vflux_gs2 = get_vflux_from_GS2(filename,tracespecs[is])

      # Interpolate fluxes to the t3core vgrid from the gs2 vgrid
      vflux_func = Spline1D(vgrid_gs2, vflux_gs2, k=1,bc="extrapolate")
      vflux[is,:] = evaluate(vflux_func,vgrid)

      df0dv_gs2 =  get_df0dv_from_GS2(filename,tracespecs[is])
      df0dv_func = Spline1D(vgrid_gs2, df0dv_gs2, k=1,bc="extrapolate")
      df0dv[is,:] = evaluate(df0dv_func,vgrid)
   end

   # Calculate Dvv and "H0" by finding linear fit of fluxes to df0/dv
   H0 = zeros(Float64,Nv)
   Dvv = zeros(Float64,Nv)
   for iv in 1:Nv
      a,b = linear_fit(df0dv[:,iv],vflux[:,iv])
      H0[iv] = a
      Dvv[iv] = -b
   end
     
   return Dvv, H0
end

function Dvv_model(vgrid)

   v0 = vmax/5.0
   D0 = 1.0

   H0 = zeros(Float64,Nv)
   Dvv = zeros(Float64,Nv)

   idx = indmin(abs(vgrid-v0))
   
   Dvv[1:idx-1] = D0

   Dvv[idx:end] = D0*(vgrid/v0).^-3

   return Dvv, H0
   
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
      push!(bulkspec, SpeciesData(mass=mass,q=charge,dens=dens,temp=temp,tprim=tprim,fprim=fprim,ispec=is) )

      # electrons
      mass = me
      charge = -qref
      dens = nref
      temp = Tref
      tprim = a / Tref
      fprim = a / nref
      push!(bulkspec, SpeciesData(mass=mass,q=charge,dens=dens,temp=temp,tprim=tprim,fprim=fprim,ispec=is) )
 
   tracespec = bulkspec[1]

   return bulkspec, tracespec
end

end

