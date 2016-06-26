module geometry
using species: filenames
using input: a, rgrid_gs2, surface_area_in, circular, rmaj, rgrid_in, grho_in, spline_k
using grids: rgrid
using NetCDF
using Dierckx

export init_geometry, surface_area, Vprime, grad_rho

Vprime = Float64[]
surface_area = Float64[]
grad_rho = Float64[]

function init_geometry()
  global Vprime, surface_area, grad_rho

  if circular
    # Use simple circular geometry

    surface_area = (rgrid*2.0*pi)*(rmaj*2.0*pi)
    Vprime = copy(surface_area)

    grho = ones(Float64,length(rgrid_in))
  else
    # Calculate geometrical properties from GS2 files.
   
    # surface_area_in is input
    if length(surface_area_in) != length(rgrid_in) 
      error("Must input surface area as a function of rgrid_in")
    end
   
    # Interpolate to find surface area and <|grad rho|> on Chebyshev grid
    area_func = Spline1D(rgrid_in, surface_area_in,k=spline_k)
    grho_func = Spline1D(rgrid_in, grho_in,k=spline_k)

    surface_area = evaluate(area_func,rgrid)
    grad_rho = evaluate(grho_func,rgrid)

    Vprime = surface_area ./ grad_rho

  end

end

end
