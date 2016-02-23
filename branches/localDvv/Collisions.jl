module Collisions
#include("Species.jl"); using .Species: SpeciesData
using Species
include("Constants.jl")

export nus, nupar

"""
nus(v::Real, a::SpeciesData, b::SpeciesData)

Returns νₛᵃᵇ(v) from Helander and Sigmar (2002), equation 3.46 (without the mass ratio factor that is cancelled in 3.40), where a and b are the test particle and target particle respectively
"""
function nus(v::Float64,a::SpeciesData,b::SpeciesData)
  local Goverv::Float64

  nuhat = nuhatvta3(a,b)
  vtb = sqrt(2.0 * b.temp / b.mass)

  if abs(v) < eps(typeof(v))
    Goverv = 2.0/(3.0*sqrt(pi))    
  else   
    Goverv = gfunc(v/vtb)/v
  end

  return nuhat* (a.mass/b.temp) * Goverv
end

"""
nupar(v::Real, a::SpeciesData, b::SpeciesData)

Returns ν∥ᵃᵇ(v) from Helander and Sigmar (2002), equation 3.47, where a and b are the test particle and target particle respectively
"""
function nupar(v::Float64,a::SpeciesData,b::SpeciesData)

  nuhat = nuhatvta3(a,b)
  vtb = sqrt(2.0 * b.temp / b.mass)

  if abs(v) < eps(typeof(v))
     return 0.0    
  else
     return nuhat* 2.0 * gfunc(v/vtb) / v^3
  end
end

function gfunc(x::Real)
  if abs(x) < eps(typeof(x))
    return zero(Real)
  else
    return (erf(x)-(2.0/sqrt(pi))*x*exp(-x^2))/(2.0*x^2)
  end
end

"""
lnlambda(a::SpeciesData, b::SpeciesData)

Calculates ln(Λᵃᵇ ) using approximations in NRL plasma formulary"
"""
function lnlambda(s1::SpeciesData, s2::SpeciesData)
  local me_threshold = 0.05*mp	# Maximum mass that can be considered an "electron", in case one needs reduced mass ratio

  # Check for electron-electron collisions
  if (s1.mass < me_threshold) && (s2.mass < me_threshold)
    Te = s1.temp/el 	
    ne = s1.dens/el
    return 23.5 - log(sqrt(ne)/Te^1.25) - sqrt(1.0e-5 + (log(Te)-2.0)^2/16.0)
  # Check for electron-ion collisions
  elseif (s1.mass < me_threshold) || (s2.mass < me_threshold)
    if (s1.mass < me_threshold) 
      me = s1.mass; Te = s1.temp/el; ne = 1.e-6*s1.dens
      mi = s2.mass; Ti = s2.temp/el; ni = 1.e-6*s2.dens; Z=s2.q/el
    else 
      mi = s1.mass; Te = s1.temp/el; ne = 1.e-6*s1.dens; Z=s1.q/el
      me = s2.mass; Ti = s2.temp/el; ni = 1.e-6*s2.dens
    end 
    if Te < Ti*Z*me/mi
      return 30.0 - log(sqrt(ni)*Z^2*(mp/mi)/Ti^1.5)
    elseif Te > 10.0*Z^2
      return 24.0 - log(sqrt(ne)/Te) 
    else
      return 23.0 - log(sqrt(ne)*Z/Te^1.5)
    end
  # Ion-ion collision
  else
    Z1 = s1.q/el; Z2 = s2.q/el
    m1 = s1.mass; m2 = s2.mass
    n1 = s1.dens; n2 = s2.dens
    T1 = s1.temp; T2 = s2.temp
    return 23.0 - log(1.e-3*el^1.5) - log( abs(Z1*Z2)* (m1+m2)*sqrt( (n1*Z1^2/T1) + (n2*Z2^2/T2)) /(m1*T2 + m2*T1 ))
  end
end

function nuhatvta3(a::SpeciesData,b::SpeciesData)
  return b.dens * b.q^2 * a.q^2 * lnlambda(a,b) / (4.0*pi*ep0^2 * a.mass^2)
end


end
