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

export wavenumbers
export compute_energy
export compute_enstrophy

abstract type AbstractSpectrumType end

struct Time <: AbstractSpectrumType end
struct Energy <: AbstractSpectrumType end
struct Helicity <: AbstractSpectrumType end

struct SpectrumData{S<:AbstractSpectrumType,
                    T, N, M<:AbstractArray{T,N}} <: AbstractArray{T,N}
    data::M
end

SpectrumData{S}(F::M) where {S<:AbstractSpectrumType,T,N,M<:AbstractArray{T,N}} = SpectrumData{S,T,N,M}(F)

Base.size(s::SpectrumData) = size(s.data)
Base.axes(s::SpectrumData) = axes(s.data)
Base.getindex(s::SpectrumData, I...) = getindex(s.data,I...)
Base.setindex!(s::SpectrumData, v, I...) = setindex!(s.data, v, I)

"""
    wavenumbers(s::SpectrumData)

Generates 1D wavenumber vector matching the size of the SpectrumData array.
"""
wavenumbers(s::SpectrumData) = collect(0:(size(s)[2] - 1))

"""
    const HydroNS

Dictionary for the HydroNS output spectra.
"""
const HydroNS = Dict(
    :time      => ("time_sp.dat",Time),
    :energy    => ("spU.dat",Energy),
    :helicity  => ("spHel.dat",Helicity),
)



include("io.jl")
include("hydroNS.jl")



end # module NSTools
