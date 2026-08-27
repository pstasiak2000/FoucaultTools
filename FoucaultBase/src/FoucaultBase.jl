module FoucaultBase

using Printf
using DelimitedFiles
using FFTW
using LinearAlgebra

export RealFFT
export ComplexFFT

export read_data
export WaveNumbers


include("fields.jl")
include("FFTs.jl")
include("io.jl")

end # module FoucaultBase
