module grids
using chebyshev: singleChebyshevDifferentiation, singleChebyshevWeights
using input: Nrad, Nv, rgrid_in, vmax, rgrid_gs2

export v,ddv,d3v,init_vgrid, init_rgrid, rgrid, ddr

v= Float64[]
ddv= Float64[]
d3v= Float64[]
rgrid = Float64[]
ddr = Float64[]

function init_rgrid()
  global rgrid
  global ddr

#  rgrid,ddr = singleChebyshevDifferentiation(Nrad, 0.0, rgrid_in[end])
#  rgrid,d3r,ddr = uniform_3pt(Nrad, 0.0, rgrid_in[end])
  rgrid,d3r,ddr = uniform_3pt(Nrad, rgrid_gs2[1], rgrid_gs2[end])
#  rgrid=rgrid[1:Nrad]
#  ddr=ddr[1:Nrad,1:Nrad]

end

function init_vgrid()
  global v, d3v, ddv

  # Set up differentiation and integration weights
#  v,d3v,ddv = uniform_3pt(Nv+1, 0.0, vmax)
  deltav = vmax/Nv
  v,d3v,ddv = uniform_3pt(Nv, 0.5*deltav, vmax+0.5*deltav)
#  v,d3v,ddv,d2dv2 = uniform_grid(Nv+1, 0.0, vmax)

  # Following scheme doesn't work because F0 varies by more than machine precision over the domain 
  # so a global derivative scheme in velocity is not appropriate
  #v,d3v = singleChebyshevWeights(Nv, 0.0, vmax)
  # Following scheme doesn't work because the fact that df/dv=0 at v=0 is needed at the second grid point also
  # Can be implemented in the future for greater accuracy
  #v,d3v,ddv,d2dv2 = uniform_grid(Nv, 0.0, vmax)

  v = v[1:Nv]
  ddv=ddv[1:Nv,1:Nv]
  d3v = (4*pi*v.^2).*d3v[1:Nv]
 
end

#> Written by Matt Landreman (scheme 12); here it's adopted to Julia
function uniform_grid(N,xmin,xmax)
  x = linspace(xmin,xmax,N)
  dx = x[2]-x[1]
  dx2 = dx^2
  w = ones(size(x)) 
  w[1] = 0.5
  w[end]=0.5
  w = w*dx
  if N<5
    error("N must be greater than 5 for this uniform-grid scheme")
  end
  D = (-diagm((1/6)*ones(N-2),2) + diagm((4/3)*ones(N-1),1) - diagm((4/3)*ones(N-1),-1) + diagm((1/6)*ones(N-2),-2))/(2*dx)
  DD = (-diagm((1/12)*ones(N-2),2) + diagm((4/3)*ones(N-1),1) - diagm((5/2)*ones(N),0) + diagm((4/3)*ones(N-1),-1) - diagm((1/12)*ones(N-2),-2))/(dx2)
  D[1,1] = -25/(12*dx)
  D[1,2] = 4/dx
  D[1,3] = -3/dx 
  D[1,4] = 4/(3*dx)
  D[1,5] = -1/(4*dx)
 
  D[2,1] = -1/(4*dx)
  D[2,2] = -5/(6*dx)
  D[2,3] = 3/(2*dx)
  D[2,4] = -1/(2*dx)
  D[2,5] = 1/(12*dx)

  D[end,end] = 25/(12*dx)
  D[end,end-1] = -4/(dx)
  D[end,end-2] = 3/(dx)
  D[end,end-3] = -4/(3*dx)
  D[end,end-4] = 1/(4*dx)

  D[end-1,end] = 1/(4*dx)
  D[end-1,end-1] = 5/(6*dx)
  D[end-1,end-2] = -3/(2*dx)
  D[end-1,end-3] = 1/(2*dx)
  D[end-1,end-4] = -1/(12*dx)

  DD[1,1] = 35/(12*dx2)
  DD[1,2] = -26/(3*dx2)
  DD[1,3] = 19/(2*dx2)
  DD[1,4] = -14/(3*dx2)
  DD[1,5] = 11/(12*dx2)

  DD[2,1] = 11/(12*dx2)
  DD[2,2] = -5/(3*dx2)
  DD[2,3] = 1/(2*dx2)
  DD[2,4] = 1/(3*dx2)
  DD[2,5] = -1/(12*dx2)

  DD[end,end] = 35/(12*dx2)
  DD[end,end-1] = -26/(3*dx2)
  DD[end,end-2] = 19/(2*dx2)
  DD[end,end-3] = -14/(3*dx2)
  DD[end,end-4] = 11/(12*dx2)

  DD[end-1,end] = 11(12*dx2)
  DD[end-1,end-1] = -5/(3*dx2)
  DD[end-1,end-2] = 1/(2*dx2)
  DD[end-1,end-3] = 1/(3*dx2)
  DD[end-1,end-4] = -1/(12*dx2)

  return x,w,D,DD

end

function uniform_3pt(N,xmin,xmax)
  x = linspace(xmin,xmax,N)
  dx = x[2] - x[1]
  w = ones(size(x)) 
  w[1] = 0.5
  w[end]=0.5
  w = w*dx

#  D = zeros(N,N)
  D = diagm(ones(N-1),1)-diagm(ones(N-1),-1)
  D[1,1] = -2.0
  D[1,2] = 2.0
  D[end,end] = 2.0
  D[end,end-1] = -2.0

  D = (1.0/(2.0*dx))*D
 
  return x,w,D

end

end
