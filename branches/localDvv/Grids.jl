module Grids
include("Input.jl")

export init_uniform_staggered_grid

function init_uniform_staggered_grid(Nv::Integer,vmax::Real)
   dv = vmax / Nv
   vgrid = collect( linspace(0.5*dv,vmax - 0.5*dv,Nv) )
   
   d3v = zeros(Nv)
   d3v[:] = dv
   return vgrid, d3v
end

end
