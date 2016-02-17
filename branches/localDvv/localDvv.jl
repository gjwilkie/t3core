include("Species.jl"); using Species
include("Constants.jl"); using Constants
include("Collisions.jl"); using Collisions
#include("DiffEq.jl"); using DiffEq
using Winston

function maxw_test()

   vts = sqrt(2.0*Tref/mass)
   vmax = sqrt(2.0*Emax/mass)
   dv = vmax / Nv

   bulkspec = SpeciesData(mass=mp,q=el,dens=1.e20,temp=T) 
   testspec = SpeciesData(mass=mp,q=el,dens=1.e20,temp=T)

   vgrid = collect(linspace(0.5*dv,vmax-dv,Nv))
   matrix = zeros(Float64,(Nv+1,Nv+1))
   source = zeros(Float64,(Nv+1))
#   matrix = zeros(Float64,(Nv,Nv))
#   source = zeros(Float64,(Nv))

   for i = 2:Nv-1
      v= vgrid[i]
      vph= 0.5*(vgrid[i]+vgrid[i+1])
      vmh= 0.5*(vgrid[i]+vgrid[i-1])
      matrix[i,i] = 0.5*( nus(vph,testspec,bulkspec)*vph^3 - nus(vmh,testspec,bulkspec)*vmh^3 - nupar(vph,testspec,bulkspec)*(vph^4/dv) - nupar(vmh,testspec,bulkspec)*(vmh^4/dv) )/dv
      matrix[i,i+1] = 0.5*( nus(vph,testspec,bulkspec)*vph^3 + nupar(vph,testspec,bulkspec)*(vph^4/dv))/dv
      matrix[i,i-1] = -0.5*( nus(vmh,testspec,bulkspec)*vmh^3 - nupar(vmh,testspec,bulkspec)*(vmh^4/dv))/dv
#println((v,nus(v,testspec,bulkspec),nupar(v,testspec,bulkspec)))
   end
   vph= dv
   vmh= 0.0
   matrix[1,1] = 0.5*( nus(vph,testspec,bulkspec)*vph^3 - nupar(vph,testspec,bulkspec)*(vph^4/dv))/dv
   matrix[1,2] = 0.5*( nus(vph,testspec,bulkspec)*vph^3 + nupar(vph,testspec,bulkspec)*(vph^4/dv))/dv
   matrix[Nv,Nv] = 1.0 * 1e20

   d3v = zeros(Float64,Nv)
   d3v[1] = 0.5*dv
   d3v[Nv] = 0.5*dv
   d3v[2:Nv-1] = dv

   d3v = 4.0*pi*d3v.*vgrid.^2

   sink = exp(-0.25*mp*vgrid.^2/T).*vgrid.^2*4.0*pi
#   sink = 0.0

   for i = 1:Nv-1
      matrix[Nv+1,i] = d3v[i]
      matrix[i,Nv+1] = sink[i]
   end    
   source[Nv+1] = 1.e20

   f = zeros(Float64,Nv+1)
   try (f = matrix\source)
   catch
      println("WARNING: Matrix is singular. Using pseudoinverse instead.")
      f = pinv(matrix)*source
   end

   fm = 1.e20*(mp/(2.0*pi*T))^1.5*exp(-vgrid.^2/vts^2)

   plot(vgrid/vts,abs(f[1:Nv]),vgrid/vts,abs(fm)  )

end

maxw_test()
