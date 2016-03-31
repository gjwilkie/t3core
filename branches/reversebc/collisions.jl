# Module that contains various routines useful for collision operators
module collisions
using constants: me,mp,el

export lnLambda, G, G_approx

# Chandrasekhar function useful for collisions against Maxwellian background
function G_approx(x::Float64)
  ep = 2.0/(3.0*sqrt(pi))
  return ep*x/(1+2*ep*x^3)
end

# Approximation to G from Stix (1975)
function G(x::Float64)
  if abs(x) < eps(typeof(x))
    return zero(Float64)
  else
    return (erf(x)-(2.0/sqrt(pi))*x*exp(-x^2))/(2.0*x^2)
  end
end

# The Coulomb logarithm, taking arguments in SI units
# Using the formula in NRL formulary
function lnLambda(; m1=mp, m2=me, Z1=1, Z2=-1, T1=1.0e3*el, T2=1.0e3*el, n1=1.0e20, n2=1.0e20)
  local me_threshold = 0.05*mp	# Maximum mass that can be considered an "electron", in case one needs reduced mass ratio

  # Check for electron-electron collisions
  if (m1 < me_threshold) && (m2 < me_threshold)
    Te = T1/el 	
    ne = n1/el
    return 23.5 - log(sqrt(ne)/Te^1.25) - sqrt(1.0e-5 + (log(Te)-2.0)^2/16.0)
  # Check for electron-ion collisions
  elseif (m1 < me_threshold) || (m2 < me_threshold)
    if (m1 < me_threshold) 
      me = m1; Te = T1/el; ne = 1.e-6*n1
      mi = m2; Ti = T2/el; ni = 1.e-6*n2; Z=Z2
    else 
      me = m2; Te = T2/el; ne = 1.e-6*n2
      mi = m1; Ti = T1/el; ni = 1.e-6*n1; Z=Z1
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
    return 23.0 - log(1.e-3*el^1.5) - log( abs(Z1*Z2)* (m1+m2)*sqrt( (n1*Z1^2/T1) + (n2*Z2^2/T2)) /(m1*T2 + m2*T1 ))
  end
end

end
