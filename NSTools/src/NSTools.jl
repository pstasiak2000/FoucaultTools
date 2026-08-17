"""
    NSTools

Tools for manipulating quantities originating from the Navier-Stokes solver.

# Overview

`NSTools` provides routines for working with scalar and vector fields, spectral outputs and global quantities.

The package is part of the `FoucaultTools` ecosystem.  
"""
module NSTools

using FoucaultBase

export read_spectrum
export read_spectrum_time

"""
    const HydroNS

Dictionary for the HydroNS output spectra.
"""
const HydroNS = Dict(
    :time => "time_sp.dat",
    :energy  => "spU.dat",
    :helicity  => "spHel.dat",
)

include("io.jl")




end # module NSTools
