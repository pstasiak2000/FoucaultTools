module FoucaultBase

using Printf
using DelimitedFiles
using FFTW


export RealFFT
export ComplexFFT

export read_data
export WaveNumbers
export dims




include("fields.jl")
include("FFTs.jl")
include("io.jl")

end # module FoucaultBase
