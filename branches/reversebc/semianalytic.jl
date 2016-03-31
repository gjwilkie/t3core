module semianalytic
using matrix: global_matrix, gindex, analytic_sd
using constants: valpha, mp, me, el, ep0
using collisions: lnLambda
using input: ir_sample, DTmix, m_trace, Z_trace, vmax
using geometry: Vprime
using grids: Nv, Nrad, v, rgrid
using species: ne, Ti, Te
using Dierckx
using PyPlot

export solve_semianalytic, F1

F1 = Float64[]

function solve_semianalytic()

  ir_set = Int64
#  plt[:style][:use]("myfig")

  # Read density from previous output file
  nalpha = zeros(Float64,Nrad)
  data = readdlm("density.dat",' ')
  rgrid_d = data[:,1]
  dens_d = data[:,2]
  dens_func = Spline1D(rgrid_d,dens_d,k=1,bc="zero")
  nhat = evaluate(dens_func,rgrid)
  nalpha = nhat.*ne

  # Find local SD at all radii
  f0sd = Array(Float64,(Nrad*Nv))
  f0sd_local = Array(Float64,(Nrad,Nv))
  for ir in 1:Nrad
    f0sd_local[ir,:] = analytic_sd(ir,nalpha[ir])
    for iv in 1:Nv
       idx = gindex(ir,iv)
       f0sd[idx] = f0sd_local[ir,iv]
    end
println(vec(f0sd_local[ir,:]))
  end

  # Calculate transport operator effect on fsd
  Tf0 = global_matrix*f0sd

  # Define interpolants
  Tf0_splines = Array(Spline1D,Nrad)
  for ir in 1:Nrad
    Tf0_local = zeros(Float64,Nv)
    for iv in 1:Nrad
      Tf0_local[iv] = Tf0[gindex(ir,iv)]
    end
    Tf0_splines[ir] = Spline1D(v,Tf0_local,k=1)
  end

  # Define 1/nupar*v^4
  function nupar_v4_inv(v_in)
    T_coll = [Ti[ir_set], Ti[ir_set], Te[ir_set]]
    n_coll = [DTmix*ne[ir_set], (1.0-DTmix)*ne[ir_set], ne[ir_set]]
    m_coll = [2.0*mp, 3.0*mp, me]
    Z_coll = [1.0,1.0,-1.0]
    nu_par_v4 = 0.0
    for is in 1:3
      # Coulomb logarithm assumes T_alpha = 10*T_i and n_alpha = 0.01*ne. Better ideas?
      logLambda = lnLambda(m1=m_trace, m2=m_coll[is], Z1=Z_trace, Z2=Z_coll[is], T1=10.0*Ti[ir_set], T2=T_coll[is], n1=0.01*ne[ir_set], n2=ne[ir_set])

      vts = sqrt(2.0*T_coll[is]/m_coll[is])

      G = erf(v_in/vts) - (2.0/sqrt(pi))*(v_in/vts).*exp(-(v_in/vts).^2)
      G = G./(2.0*(v_in/vts).^2)

      # Construct (nuhat*vta^3) for this species pair
      # Eq. (3.48) of Hellander and Sigmar, times vta^3
      nuhat_vta3 = n_coll[is]*Z_coll[is]^2*Z_trace^2*el^4*logLambda
      nuhat_vta3 = nuhat_vta3 / (4.0*pi*(ep0*m_trace)^2)

      nu_par_v4 = nu_par_v4 + 2.0*nuhat_vta3*G*v_in

    end   
    return 1.0/nu_par_v4
  end
 
  # Define nu_s/(nupar*v)
  function nus_nupar_v(v_in)
    T_coll = [Ti[ir_set], Ti[ir_set], Te[ir_set]]
    n_coll = [DTmix*ne[ir_set], (1.0-DTmix)*ne[ir_set], ne[ir_set]]
    m_coll = [2.0*mp, 3.0*mp, me]
    Z_coll = [1.0,1.0,-1.0]
    nu_s_v2 = 0.0
    nu_par_v3 = 0.0
    for is in 1:3
      # Coulomb logarithm assumes T_alpha = 10*T_i and n_alpha = 0.01*ne. Better ideas?
      logLambda = lnLambda(m1=m_trace, m2=m_coll[is], Z1=Z_trace, Z2=Z_coll[is], T1=10.0*Ti[ir_set], T2=T_coll[is], n1=0.01*ne[ir_set], n2=ne[ir_set])

      vts = sqrt(2.0*T_coll[is]/m_coll[is])

      G = erf(v_in/vts) - (2.0/sqrt(pi))*(v_in/vts).*exp(-(v_in/vts).^2)
      G = G./(2.0*(v_in/vts).^2)

      # Construct (nuhat*vta^3) for this species pair
      # Eq. (3.48) of Hellander and Sigmar, times vta^3
      nuhat_vta3 = n_coll[is]*Z_coll[is]^2*Z_trace^2*el^4*logLambda
      nuhat_vta3 = nuhat_vta3 / (4.0*pi*(ep0*m_trace)^2)

      nu_s_v2 = nu_s_v2 + nuhat_vta3*m_trace*(G*v_in)/T_coll[is]

      nu_par_v3 = nu_par_v3 + 2.0*nuhat_vta3*G

    end
    return nu_s_v2/nu_par_v3
  end

function inner_integrand(v_in)
  return evaluate(Tf0_splines[ir_set],v_in)/Vprime[ir_set]
end

function outer_integrand(v_in)
  integral,err = quadgk(inner_integrand,0.0,v_in)
  return 2.0*exp(2.0*nus_nupar_v(v_in))*nupar_v4_inv(v_in)*integral
end

  ir_set = 1
  test = inner_integrand(v)
  plot(v/valpha,test)
  savefig("test.png")
  cla()
  clf()

  F1 = zeros(Float64,(Nrad,Nv))
  for ir in 1:Nrad
    ir_set = ir
    for iv in 1:Nv
      integral,err = quadgk(outer_integrand,v[iv],vmax,reltol=1.e-4)
      F1[ir,iv] = exp(-2.0*nus_nupar_v(v[iv]))*integral
    end
  end

  semilogy(v/valpha,abs(vec(f0sd_local[ir_sample,:])))
  semilogy(v/valpha,abs(vec(f0sd_local[ir_sample,:]+F1[ir_sample,:])))
#  plot(v/valpha,abs(vec(f0sd_local[ir_sample,:])))
#  plot(v/valpha,abs(vec(f0sd_local[ir_sample,:]+F1[ir_sample,:])))
  legend(("SD","with perturb."))
  savefig("F1.png")
  cla()
  clf()

end

end
