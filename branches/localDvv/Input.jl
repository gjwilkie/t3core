include("Constants.jl"); using Constants

# Simulation parameters
const Nv = 400             # Resolution in energy
const tavg_frac = 0.5      # Fraction of runtime over which to time-average

# Normalization constants with respect to GS2 (SI units)
const dens = 1.e20
const mref = mp
const qref = el
const Tref = 10.0e3*el
const vmax = 6.0*sqrt(2.0*Tref/mref)
const a = 1.0
const rhostar = 0.001

const filenames = ["slowmode.out.nc"]
const traceidx = [1,3]
const bulkidx = [1,2]

const agk = true

# Derived quantities
const vref = sqrt(2.0*Tref/mref)
