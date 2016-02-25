include("Constants.jl")

# Simulation parameters
const Nv = 800             # Resolution in energy
const tavg_frac = 0.5      # Fraction of runtime over which to time-average

# Boundary condition options
const nedge = 1.0e20

# Normalization constants with respect to GS2 (SI units)
const nref = 1.e20
const mref = mp
const qref = el
const Tref = 10.0e3*el
const vmax = 8.0*sqrt(2.0*Tref/mref)
const a = 1.0
const rhostar = 0.001

const filenames = ["slowmode.out.nc"]
const traceidx = [1,3]
const bulkidx = [1,2]

const turb_rescale = 0.0

const agk = true

# Derived quantities
const vref = sqrt(2.0*Tref/mref)
