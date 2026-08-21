### ====================================================== ##
### FFTs (variables and plans)
### ====================================================== ##

abstract type FFTKind end

struct RealFFT <: FFTKind end
struct ComplexFFT <: FFTKind end

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

### By default the WaveNumbers relate to real fields

Base.size(kk::WaveNumbers) = kk.N
dims(kk::WaveNumbers{D}) where D = D
fft_kind(kk::WaveNumbers{D,T,F}) where {D,T,F<:FFTKind} = F()

