#!/bin/julia
using Collisions
using Grids
using Input
using Species
using Winston

function main()

   info("Preparing data...")
   vgrid, d3v = init_uniform_staggered_grid(Nv,vmax)

   dv = vgrid[2] -vgrid[1]
   vgrid_ph = vgrid + 0.5*dv
   vgrid_mh = vgrid - 0.5*dv

#   bulkspecs,tracespecs= get_species_from_GS2(filenames[1])
   bulkspecs,tracespecs= get_test_species()

#   H0, Dvv = calculate_Dvv(filenames[1],tracespecs,vgrid)
#   H0_ph, Dvv_ph = calculate_Dvv(filenames[1],tracespecs,vgrid_ph)
#   H0_mh, Dvv_mh = calculate_Dvv(filenames[1],tracespecs,vgrid_mh)
   H0, Dvv = Dvv_model(vgrid)
   H0_ph, Dvv_ph = Dvv_model(vgrid_ph)
   H0_mh, Dvv_mh = Dvv_model(vgrid_mh)

   matrix = zeros(Float64,(Nv+1,Nv+1))
   source = zeros(Float64,(Nv+1))

   info("Populating matrix...")
   for i = 2:Nv-1
      v= vgrid[i]
      vph= 0.5*(vgrid[i]+vgrid[i+1])
      vmh= 0.5*(vgrid[i]+vgrid[i-1])

      # Collision operator
      for is in length(bulkspecs)
         matrix[i,i] += 0.5*( nus(vph,tracespecs[1],bulkspecs[is])*vph^3 - nus(vmh,tracespecs[1],bulkspecs[is])*vmh^3 - nupar(vph,tracespecs[1],bulkspecs[is])*(vph^4/dv) - nupar(vmh,tracespecs[1],bulkspecs[is])*(vmh^4/dv) )/dv
         matrix[i,i+1] += 0.5*( nus(vph,tracespecs[1],bulkspecs[is])*vph^3 + nupar(vph,tracespecs[1],bulkspecs[is])*(vph^4/dv))/dv
         matrix[i,i-1] += -0.5*( nus(vmh,tracespecs[1],bulkspecs[is])*vmh^3 - nupar(vmh,tracespecs[1],bulkspecs[is])*(vmh^4/dv))/dv
      end

      # Turbulence
      matrix[i,i] += 0.5*( -H0_ph[i]*vph^2 + H0_mh[i]*vmh^2 - 2.0*Dvv_ph[i]*(vph^2/dv) - 2.0*Dvv_mh[i]*(vmh^2/dv) )/dv
      matrix[i,i+1] += 0.5*( -H0_ph[i]*vph^2 + 2.0*Dvv_ph[i]*(vph^2/dv))/dv
      matrix[i,i-1] += 0.5*( H0_mh[i]*vmh^2 + 2.0*Dvv_mh[i]*(vmh^2/dv))/dv

   end

   info("Setting boundary conditions...")
   # Boundary cases
   vph= dv
   vmh= 0.0

   for is in length(bulkspecs)
      matrix[1,1] += 0.5*( nus(vph,tracespecs[1],bulkspecs[is])*vph^3 - nupar(vph,tracespecs[1],bulkspecs[is])*(vph^4/dv))/dv
      matrix[1,2] += 0.5*( nus(vph,tracespecs[1],bulkspecs[is])*vph^3 + nupar(vph,tracespecs[1],bulkspecs[is])*(vph^4/dv))/dv
   end

   matrix[1,1] += 0.5*( -H0_ph[1]*vph^2 - 2.0*Dvv_ph[1]*(vph^2/dv))/dv
   matrix[1,2] += 0.5*( -H0_ph[1]*vph^2 + 2.0*Dvv_ph[1]*(vph^2/dv))/dv

   if false
      matrix[Nv,Nv] = nedge
   else
      vph= vmax
      vmh= vmax-dv

      for is in length(bulkspecs)
         matrix[Nv,Nv] += 0.5*( - nus(vmh,tracespecs[1],bulkspecs[is])*vmh^3 - nupar(vmh,tracespecs[1],bulkspecs[is])*(vmh^4/dv) )/dv
         matrix[Nv,Nv-1] += -0.5*( nus(vmh,tracespecs[1],bulkspecs[is])*vmh^3 - nupar(vmh,tracespecs[1],bulkspecs[is])*(vmh^4/dv))/dv
      end

      matrix[Nv,Nv] += 0.5*( H0_mh[Nv]*vmh^2 - 2.0*Dvv_mh[Nv]*(vmh^2/dv) )/dv
      matrix[Nv,Nv-1] += 0.5*( H0_mh[Nv]*vmh^2 + 2.0*Dvv_mh[Nv]*(vmh^2/dv))/dv
   end

   sink = exp(-0.25*mref*vgrid.^2/Tref).*vgrid.^2*4.0*pi
#   sink = 0.0

   for i = 1:Nv-1
      matrix[Nv+1,i] = d3v[i]
      matrix[i,Nv+1] = sink[i]
   end    
   source[Nv+1] = nedge

   info("Inverting matrix...")
   f = zeros(Float64,Nv+1)
   try (f = matrix\source)
   catch
      println("WARNING: Matrix is singular. Using pseudoinverse instead.")
      f = pinv(matrix)*source
   end

   vts = sqrt(2.0*Tref/mref)
   fm = nedge*(mref/(2.0*pi*Tref))^1.5*exp(-vgrid.^2/vts^2)

#   plot(vgrid/vts,abs(f[1:Nv]),vgrid/vts,abs(fm)  )
   semilogy(vgrid/vts,abs(f[1:Nv]),vgrid/vts,abs(fm)  )


end

function maxw_test()

   vts = sqrt(2.0*Tref/mass)
   vmax = sqrt(2.0*Emax/mass)
   dv = vmax / Nv

   bulkspec = SpeciesData(mass=mp,q=el,dens=1.e20,temp=T) 
   testspec = SpeciesData(mass=mp,q=el,dens=1.e20,temp=T)

   vgrid = collect(linspace(0.5*dv,vmax-dv,Nv))
   matrix = zeros(Float64,(Nv+1,Nv+1))
   source = zeros(Float64,(Nv+1))

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
   matrix[Nv,Nv] = nedge

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
   source[Nv+1] = nedge

   f = zeros(Float64,Nv+1)
   try (f = matrix\source)
   catch
      println("WARNING: Matrix is singular. Using pseudoinverse instead.")
      f = pinv(matrix)*source
   end

   vts = sqrt(2.0*Tref/mref)
   fm = 1.e20*(mp/(2.0*pi*T))^1.5*exp(-vgrid.^2/vts^2)

   plot(vgrid/vts,abs(f[1:Nv]),vgrid/vts,abs(fm)  )

end

main()
