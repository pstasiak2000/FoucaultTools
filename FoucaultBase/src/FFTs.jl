### ====================================================== ##
### FFTs (variables and plans)
### ====================================================== ##

abstract type FFTKind end

"""
    RealFFT

Specifies FFT to use real-to-complex ffts.
"""
struct RealFFT <: FFTKind end

"""
    ComplexFFT

Specifies FFT to use complex-to-complex ffts.
"""
struct ComplexFFT <: FFTKind end


"""
    WaveNumbers{D,T,F}(N, L, x)

Store the wavenumber information associated with a `D`-dimensional field.

# Type Parameters

- `D`: Number of spatial dimensions.
- `T`: Floating-point element type used for the domain lengths and wavenumbers.
- `F <: FFTKind`: FFT representation used to construct the wavenumbers, e.g.
  [`RealFFT`](@ref) or [`ComplexFFT`](@ref).

# Fields

- `N::NTuple{D,Int}`: Number of grid points in each spatial dimension.
- `L::NTuple{D,T}`: Physical size of the domain in each spatial dimension.
- `x::NTuple{D,Vector{T}}`: Wavenumber vectors for each spatial dimension.

The number of entries in each wavenumber vector may differ from the
corresponding value in `N` when a real FFT is used, since only the
non-negative frequencies are retained in the first dimension.
"""
struct WaveNumbers{D,T,F<:FFTKind}
    N::NTuple{D,Int}
    L::NTuple{D,T}
    x::NTuple{D,Vector{T}}
end

"""
    WaveNumbers{ComplexFFT}(N::NTuple{D,Int}, L::NTuple{D,T}) where {D,T}

Create wavenumbers based on a `D`-dimensional *complex field* of size `N`, characterised by a box size with side length `L`.

Wavenumbers are determined by [AbstractFFTs.jl](https://juliamath.github.io/AbstractFFTs.jl/stable/api/#FFT-and-FFT-planning-functions) using `fftfreq`.
Further information on how the FFTs are constructed can be found [here](https://www.matecdev.com/posts/fft-tutorial-basics.html).
"""
function WaveNumbers{ComplexFFT}(N::NTuple{D,Int}, L::NTuple{D,T}) where {D,T}
    x = ntuple(i -> fftfreq(N[i], 2π * N[i] / L[i]), Val(D))
    Nk = ntuple(i -> length(x[i]), Val(D))
    return WaveNumbers{D,T,ComplexFFT}(Nk,L,x)
end


"""
    WaveNumbers{RealFFT}(N::Tuple{D,Int}, L::NTuple{D,T})

Construct wavenumbers for a `D`-dimensional field represented using a
real-to-complex FFT.

The first dimension uses `rfftfreq`.
"""
function WaveNumbers{RealFFT}(N::NTuple{D,Int}, L::NTuple{D,T}) where {D,T} 
        x = (rfftfreq(N[1], 2π * N[1] / L[1]),
             ntuple(i -> fftfreq(N[i+1], 2π * N[i+1] / L[i+1]), Val(D - 1))...,
            )
    Nk = ntuple(i -> length(x[i]), Val(D))
    return WaveNumbers{D,T,RealFFT}(Nk,L,x)
end

WaveNumbers{F}(N::NTuple{D,Int}) where {D,F<:FFTKind} = WaveNumbers{F}(N, ntuple(_ -> 2π, Val(D)))
WaveNumbers(N::NTuple{D,Int}) where D = WaveNumbers{RealFFT}(N)

"""
    WaveNumbers{F}(N::Vararg{Int}) where {F<:FFTKind}

Specify the wavenumbers based on the dimensionality set by the number of input parameters `length(N)`.

The default behaviour is to use `RealFFT` as the FFT type.
"""
WaveNumbers{F}(N::Vararg{Int}) where {F<:FFTKind} = WaveNumbers{F}(N)
WaveNumbers(N::Vararg{Int}) = WaveNumbers{RealFFT}(N)


Base.size(kk::WaveNumbers) = kk.N
dims(kk::WaveNumbers{D}) where D = D
fft_kind(kk::WaveNumbers{D,T,F}) where {D,T,F<:FFTKind} = F()

