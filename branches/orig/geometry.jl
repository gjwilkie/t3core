module geometry
using species: filenames
using input: a, rgrid_gs2, surface_area_in, circular, rmaj, rgrid_in, grho_in, spline_k
using grids: rgrid
using NetCDF
using Dierckx

export init_geometry, surface_area, Vprime, Vprime_global, surface_area_global, grho_global

Vprime = Float64[]
Vprime_global = Float64[]
surface_area = Float64[]
surface_area_global = Float64[]
grho_global = Float64[]

function init_geometry()
  global Vprime, surface_area, Vprime_global, surface_area_global, grho_global

  if circular
    # Use simple circular geometry

    surface_area = (rgrid*2.0*pi)*(rmaj*2.0*pi)
    Vprime = copy(surface_area)

    surface_area_global = (rgrid_in*2.0*pi)*(rmaj*2.0*pi)
    Vprime_global = copy(surface_area_global)
    grho_global = ones(Float64,length(rgrid_in))
  else
    # Calculate geometrical properties from GS2 files.
   
    # surface_area_in is input
    if length(surface_area_in) != length(rgrid_in) 
      error("Must input surface area as a function of rgrid_in")
    end
   
    # grho input now, no longer computed because it is needed where gs2 wasn't ran
#    grad_rho_gs2 = Float64[]
#    for file in filenames
#      grho = ncread(file,"grho")[:] 
#      theta = ncread(file,"theta")[:] 
#      bmag = ncread(file,"bmag")[:] 
#      gradpar = ncread(file,"gradpar")[:] 
#      Ntheta=length(theta)
#      delthet = theta[2:Ntheta] - theta[1:Ntheta-1]
#      push!(delthet,0.0) 
#
#      # Save the flux-surface average of |grad(rho)| to grad_rho_gs2
#      push!(grad_rho_gs2,sum(grho.*delthet./(a*bmag.*gradpar)))
#    end 

    # Interpolate to find surface area and <|grad rho|> on Chebyshev grid
    area_func = Spline1D(rgrid_in, surface_area_in,k=spline_k)
    grho_func = Spline1D(rgrid_in, grho_in,k=spline_k)

    surface_area = evaluate(area_func,rgrid)
    grad_rho = evaluate(grho_func,rgrid)

    surface_area_global = copy(surface_area_in)
    grho_global = copy(grho_in)

    Vprime = surface_area ./ grad_rho

#    println("Vprime = ",Vprime)
#    println("(2*pi*r)^2*R = ",(2.0*pi*rgrid).^2)*rmaj

    Vprime_global = surface_area_global ./ grho_in
#    Vprime_global = copy(surface_area_global) 
  end

end

end
