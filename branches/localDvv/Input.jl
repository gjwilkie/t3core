include("Constants.jl")

# Simulation parameters
const Nv = 800             # Resolution in energy
const tavg_frac = 0.1      # Fraction of runtime over which to time-average

# Normalization constants with respect to GS2 (SI units)
const nref = 1.e6
const mref = mp
const qref = el
const Tref = 1.0e1*el
const rhostar = 0.1
const a = 1.0e7
const ntot = 1.0e4
#const nref = 1.e20
#const mref = 2.0*mp
#const qref = el
#const Tref = 10.0e3*el
#const rhostar = 0.001
#const a = 2.0
#const ntot = 1.0e18

const vmax = 6.0*sqrt(2.0*Tref/mref)

#const filenames = ["r1.nc"]
#const traceidx = [4,5]
#const bulkidx = [1,2,3]
const filenames = ["oxygen.out.nc"]
const traceidx = [3,4]
const bulkidx = [1,2]

const turb_rescale = 1.0

const agk = true

# Derived quantities
const vref = sqrt(2.0*Tref/mref)
