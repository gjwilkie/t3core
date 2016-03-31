module input
using constants
using Dierckx

export read_input,Nrad,Nv,vmax,Nrad_gs2,tavg,deltat,tracespecs,mref,a,rhostar, rgrid_gs2, nedge, Tashfac, qref, Nt, rgrid_in, Te_in, Ti_in, ne_in, Z_trace, m_trace, rmaj, diffmodel, D0, ir_sample, dilution_model, brv, bvr, vflux_fac, semianalytic_on, zerosource, ejection_mode, diff_power, diff_v0, diff_D0, change_diffmodel, maxwellian_edge, surface_area_in, grho_in, constantD,turbfac, emrescale, spline_k, dilute_fac, ashmode, vt_temp_fac

Nrad=Int64
Nrad_gs2=Int64
Nv=Int64
nedge=Float64
tavg=Float64[]
Tashfac=Float64
deltat=Float64
Nt=Int64
vmax=Float64
ne_in=Float64
Te_in=Float64
Ti_in=Float64
rgrid_in=Float64
mref=Float64
qref=Float64
m_trace=Float64
Z_trace=Float64
a=Float64
rmaj=Float64
rhostar=Float64
circular=Bool
maxwellian_edge = Bool
semianalytic_on=Bool
diffmodel=Int64
ir_sample=Int64
dilution_model=Int64
DTmix=Float64
tracespecs=Int64[]
rgrid_gs2=Float64[]
D0 = Float64[]
brv = Float64[]
bvr = Float64[]
vflux_fac = Float64[]
zerosource = Bool
emrescale = Bool
ejection_mode = Bool
diff_D0 = Float64[]
diff_v0 = Float64[]
diff_power = Float64[]
surface_area_in = Float64[]
grho_in = Float64[]
constantD = Float64[]
turbfac = Float64[]
spline_k = Int64
dilute_fac = Float64[]
vt_temp_fac = Float64[]
ashmode = Bool

function read_input()
  global nedge, Nv, Nrad, circular, Tashfac, deltat, tracespecs,vmax,mref,qref,a,rhostar, rgrid_gs2, tavg, deltat, Nrad_gs2, Nt, rgrid_in, Te_in, Ti_in, ne_in, DTmix, m_trace, Z_trace, rmaj, diffmodel, ir_sample, dilution_model, vflux_fac, semianalytic_on, ash_cutoff, ash_accuracy, zerosource, ejection_mode, diff_power, diff_v0, diff_D0, maxwellian_edge, surface_area_in, grho_in, constantD, turbfac, emrescale, spline_k, dilute_fac, ashmode, vt_temp_fac
  turbfac=1.0e-3
  nedge=1.e19         		 # Edge density (in m^-3)
  maxwellian_edge = false

# Resolution and domain:
  Nv=400			# Number of speed grid points
  vmax = 1.05*sqrt(2.0*Ealpha/(4.0*mp))
#  vmax = 1.05*sqrt(2.0*Ealpha/(4.0*mp))

  diffmodel=1
  # 1 = All four diffusion coefficients
  # 2 = Only radial transport
  # 3 = Only radial diffusion 
  # -1= Scaling based on results from Wilkie,JPP 2015
  # -2= Constant D throughout entire domain
  constantD=0.1
  vflux_fac = 1.0
#  diff_v0 = 1.5			# Multiple of Helium thermal speed (at Ti) at which scaling with energy starts
#  diff_power = -3.0			# Power by which diffusion coefficient scales with speed
#  diff_D0 = 0.2			# Constant diffusion at low energy, as multiple of rhostar^2*vti*a
  diff_v0 = 1.5			# Multiple of Helium thermal speed (at Ti) at which scaling with energy starts
  diff_power = -1.0			# Power by which diffusion coefficient scales with speed
  diff_D0 = 0.02			# Constant diffusion at low energy, as multiple of rhostar^2*vti*a

#  vflux_fac = 1.0e-2

  dilution_model = 0
  # 0 = No dilution effec
  # 1 = Read from density.dat if it exists and adjust turbulent amplitude according to ITG scaling. Write same file when complete.
  # 2 = Use same density.dat file, but adjust fusion source instead
  # 3 = Dilution effect on both source and turbulent amplitude
  dilute_fac = 6.0

  ejection_mode=false

  ashmode = false
  
  semianalytic_on = false

  zerosource = false

  circular=false		# Use circular flux surfaces, regardless of what GS2 says
 
  emrescale=false

  Nrad=30  		# Radial grid resolution. Not necessarily the same as the number of GS2 simluations.

  Tashfac=1.0            # Multiplicative factor to determine edge ash temperature from average species temperature

  deltat= -0.5             # Timestep in s for non-steady-state solution. Set as negative for steady-state
  Nt = 10

  # Parameters used to make sense of GS2 data:

  tavg = 200.0		# Time over which to time-average. All time-averages will be the last tavg units (in a/vref)

  tracespecs=[4,5]	# Indices of the species from which to construct a generalized flux. Must both have same mass and charge in GS2 sims

  # Reference quantities
  # **** Assumed to be the same for all GS2 files! (this is contrary to what Trinity does) ****
  mref = 2.0*mp  		# Reference mass (in kg)
  qref = el		# Reference charge (in Coulombs, usually the elementary charge)*vref
  a = 2.0                 # Minor radius in meters. Also the reference length of GS2 sims
  rmaj = 6.0		# Major radius, required for cylindrical geometry
  DTmix = 0.5		# Fraction of ions that are Deuterium (rest are Tritium)

  m_trace = 4.0*mp
  Z_trace = 2.0
#  m_trace = 2.0*mp
#  Z_trace = 1.0

  vt_temp_fac = 1.0
 
  ir_sample = 1
  
  # Grid that determines background profiles (which are used to calculate alpha profile regardless if GS2 is run there or not)
  rgrid_in = [0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]*a
  ne_in = 1.0e20*[1.009,1.009,1.009,1.009,1.009,1.008,1.008,1.007,1.006]
  Te_in = 1000.0*el*[23.49,23.18,22.26,20.73,18.60,15.95,12.93,9.74,6.63]
  Ti_in = 1000.0*el*[19.49,19.24,18.49,17.26,15.54,13.39,10.94,8.36,5.84]
#  Te_in = 1000.0*el*10.0*ones(Float64,9)
#  Ti_in = 1000.0*el*10.0*ones(Float64,9)
  surface_area_in = [0.0,60.217,120.59,182.3,245.67,310.38,375.43,440.67,507.55]
  grho_in = (1.0/a)*[.8671,.869,.8591,.8390,.8209,.8128,.8088,.8,.7682]
  
#  tavg = [1500.0,800.0,300.0,150.0]
# For ash sims:
  tavg = [1500.0,500.0,200.0,110.0]

# For comparison with Sigmar, Gormley and Kamelander:
#  Te_in = Te_in*(10000.0*el/Te_in[6])
#  Ti_in = Te_in
#  ne_in = 2.0*ne_in

  # The radii (in m) at which the GS2 runs are made
  rgrid_gs2 = a*[0.5,0.6,0.7,0.8]
  rhostar = [0.00208,0.00184,0.00158,0.00130] 
#o  rhostar = rhostar.*sqrt([0.333,0.666,2.0,3.0])
#i  rhostar = rhostar.*sqrt([3.0,2.0,0.666,0.333])
 
  spline_k=3 			# Order of spline to use for radial dependence (velocity dependence is always linear)

  # Surface area (in m^2) of the flux surface at each radius
  # If not defined, will be automatically calcualted from cylindrical geometry

  Nrad_gs2 = length(rgrid_gs2)
  
#  if !isdefined(:surface_area_gs2)
#    circular = true
#  end

  if (length(rhostar)) < length(rgrid_gs2)
    rhostar = fill(rhostar[1],length(rgrid_gs2))
  end

  if (length(tavg)) < length(rgrid_gs2)
    tavg = fill(tavg[1],length(rgrid_gs2))
  end

#  if !isdefined(:ir_sample)
#    ir_sample = int(Nrad/2)
#  end

end

function input_for_analytic_test()
  global a, circular, vmax, rgrid_gs2, Nrad, Nv, rgrid_in, deltat, brv,bvr, D0, rmaj, semianalytic_on, surface_area_in, spline_k
  brv = 0.5
  bvr = 0.5
#  b = 0.0
  D0 = 5.0
  a = 1.0
  rmaj = 3.0
  circular = true
  vmax = 5.0*valpha
  rgrid_gs2 = a*[0.3,1.0]
  Nrad = 10
  Nv = 100
#  rgrid_in = [0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]*a
  rgrid_in = [0.0:0.01:1.0]*a
  deltat = -1.0
  spline_k = 3
  semianalytic_on = false
end

function change_diffmodel(diff_in)
  global diffmodel
  diffmodel = diff_in
end

end

#----------------------------

