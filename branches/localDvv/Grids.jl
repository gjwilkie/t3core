module Grids
include("Input.jl")

type Grid
   pts::Array{Float64,1}
   intwgt::Array{Float64,1}
end
   
function init_uniform_staggered_grid(Nv::Integer,vmax::Real)
   
   return vgrid
end
