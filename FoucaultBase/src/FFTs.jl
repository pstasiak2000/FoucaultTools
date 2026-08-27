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
    WaveNumbers{D,T,F}

Store the wavenumber information associated with a `D`-dimensional field.

# Type Parameters

- `D`: Number of spatial dimensions.
- `T`: Floating-point element type used for the domain lengths and wavenumbers.
- `F <: FFTKind`: FFT representation used to construct the wavenumbers, e.g.
  [`RealFFT`](@ref) or [`ComplexFFT`](@ref).

# Fields

- `N::NTuple{D,Int}`: Size of the wavenumber tuple.
- `L::NTuple{D,T}`: Physical size of the domain in each spatial dimension.
- `x::NTuple{D,Vector{T}}`: Wavenumber vectors for each spatial dimension.

# Construction

When constructing the `WaveNumbers` structure as below: 
 
```julia
WaveNumbers{D,T,F}(Nx, L)
WaveNumbers{F<:FFTKind)(Nx,L)
WaveNumbers(Nx::Vararg{Int})
WaveNumbers{Nx::Tuple{D,Int})
```

it is important to note that the argument `Nx` is related to the size of the physical field related to the wavenumbers which you are trying to create.
For example, if a field `f` has size `(128,128,128)`, then you must pass `Nx=(128,128,128)` as the relevant parameter to construct. 
If this field is real, then the associated wavenumber space given by the structure field `N` should have size `(65,128,128)`.
See [`size`](@ref) for more details.

# Notes 

The FFT representation dictates the construction of the wavenumbers for the first element in the wavenumbers tuple.
When selecting [`RealFFT`](@ref) the method `rfftfreq` is used instead of the general `fftreq` to compute the first wavenumbers.

The number of entries in each wavenumber vector may differ from the
corresponding value in `N` when a real FFT is used, since only the
non-negative frequencies are retained in the first dimension.

Wavenumbers are determined by [AbstractFFTs.jl](https://juliamath.github.io/AbstractFFTs.jl/stable/api/#FFT-and-FFT-planning-functions) using `fftfreq`.
Further information on how the FFTs are constructed can be found [here](https://www.matecdev.com/posts/fft-tutorial-basics.html).
"""
struct WaveNumbers{D,T,F<:FFTKind}
    N::NTuple{D,Int}
    L::NTuple{D,T}
    x::NTuple{D,Vector{T}}
end

function WaveNumbers{ComplexFFT}(N::NTuple{D,Int}, L::NTuple{D,T}) where {D,T}
    x = ntuple(i -> fftfreq(N[i], 2π * N[i] / L[i]), Val(D))
    Nk = ntuple(i -> length(x[i]), Val(D))
    return WaveNumbers{D,T,ComplexFFT}(Nk,L,x)
end

function WaveNumbers{RealFFT}(N::NTuple{D,Int}, L::NTuple{D,T}) where {D,T} 
        x = (rfftfreq(N[1], 2π * N[1] / L[1]),
             ntuple(i -> fftfreq(N[i+1], 2π * N[i+1] / L[i+1]), Val(D - 1))...,
            )
    Nk = ntuple(i -> length(x[i]), Val(D))
    return WaveNumbers{D,T,RealFFT}(Nk,L,x)
end

WaveNumbers{F}(N::NTuple{D,Int}) where {D,F<:FFTKind} = WaveNumbers{F}(N, ntuple(_ -> 2π, Val(D)))
WaveNumbers(N::NTuple{D,Int}) where D = WaveNumbers{RealFFT}(N)

WaveNumbers{F}(N::Vararg{Int}) where {F<:FFTKind} = WaveNumbers{F}(N)
WaveNumbers(N::Vararg{Int}) = WaveNumbers{RealFFT}(N)


"""
    Base.size(w::WaveNumbers)

Return the size of the Fourier-space grid represented by `w`.

For a complex FFT, this is the same as the physical grid size. For a
real FFT, the first dimension is reduced to `N[1] ÷ 2 + 1`, consistent
with the output of [`rfft`](@ref).

# Examples

```julia
w = WaveNumbers{RealFFT}(128, 128, 128)
```

size(w)
# (65, 128, 128)
"""
Base.size(w::WaveNumbers) = kk.N
Base.ndims(w::WaveNumbers{D}) where D = D
fft_kind(w::WaveNumbers{D,T,F}) where {D,T,F<:FFTKind} = F()

