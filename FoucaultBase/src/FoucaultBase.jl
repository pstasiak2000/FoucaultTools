module FoucaultBase

using Printf
using DelimitedFiles
using FFTW

include("FFTs.jl")
include("io.jl")

export RealFFT
export ComplexFFT

export WaveNumbers
export dims


end # module FoucaultBase
