include("Constants.jl"); using Constants

# Simulation parameters
const Nv = 400             # Resolution in energy
const tavg_frac = 0.5      # Fraction of runtime over which to time-average

# Normalization constants with respect to GS2 (SI units)
const dens = 1.e20
const mref = mp
const qref = el
const Tref = 10.0e3*el
const Emax = 18.0*Tref
const a = 1.0

const filenames = ["slowmode.out.nc"]
const traceidx = [1,3]
const bulkidx = [1,2]
