### ====================================================== ##
### HydroNS Routines
### ====================================================== ##
# 
#  This file contains the routines to be used for analysis and processing of the spectrum data primarily found in the HydroNS/ directory in the outout of FOUCAULT
#
#
#

"""
    compute_energy(s::SpectrumData{NSTools.Energy})

Computes the kinetic energy ```E_{kin}``` evolution from the energy spectrum data.

# Formula
```math
    E_{kin} = \\int_0^{\\infty} E(k) dk
```
"""
compute_energy(s::SpectrumData{Energy}) = vec(sum(s,dims=2))

"""
    compute_enstrophy(s::SpectrumData{NSTools.Energy})

Computes the enstrophy ```\\omega^2``` evolution from the energy spectrum data. 

# Formula
```math
    \\omega^2} = \\int_0^{\\infty} k^2 E(k) dk
```

where ```\\nu``` is the numerical viscosity of the fluid.
"""
function compute_enstrophy(s::SpectrumData{Energy}) 
    kk = wavenumbers(s)
    return 2.0 .* (s * (kk.^2)) 
end

"""
    compute_turb_KE(s::SpectrumData{Energy})

Computes the turbulent kinetic energy (energy of the fluctuating fluid per unit mass) from the energy spectrum data.
"""
@inline compute_turb_KE(s::SpectrumData{Energy}) = sum(view(s,:,2:size(s)[2]),dims=2) 


"""
    compute_isotropic_rms(s::SpectrumData{Energy})

Computes the isotropic root mean square of the velocity using [`compute_turb_KE`](@ref).

# Formula
```math
u_{rms} = \\sqrt{\\langle u'^2 \\rangle} = \\sqrt{(2/3) k'}
```
"""
@inline compute_isotropic_rms(s::SpectrumData{Energy}) = sqrt.((2.0/3.0) .* compute_turb_KE(s))

"""
    const ApproxRMS

Dictionary for the methods to approximate the root mean square velocity.
"""
const ApproxRMS = Dict(  
    :isotropic => s -> compute_isotropic_rms(s),
)

"""
    compute_rms_velocity(s::SpectrumData{Energy} ; rms_type::Symbol)

Computes the root mean square velocity from the energy spectrum data. 

By default it uses the isotropic approximation [`compute_isotropic_rms`](@ref) to compute the turbulent kinetic energy. 
"""
compute_rms_velocity(s::SpectrumData{Energy} ; rms_type::Symbol=:isotropic) = ApproxRMS[rms_type](s)


"""
    compute_integral_scale(s::SpectrumData{Energy} ; rms_type::Symbol)

Computes the integral length scale ```L_0``` from the energy spectrum data.

The default behaviour of `rms_type` is to use the isotropic approximation [`compute_isotropic_rms`](@ref) to compute the turbulent kinetic energy.

# Formula
```math
L_0 = \\frac{\\pi}{2U_{rms}^2} \\int_{k_0}^{\\infty} \\frac{E(k)}{k} dk
```
where ```k_0``` is the smallest wavenumber greater than 0. 
"""
function compute_integral_scale(s::SpectrumData{Energy} ; kwargs...)
    kk = wavenumbers(s)
    u_rms_sq = compute_rms_velocity(s; kwargs...) .^ 2 .+ eps(eltype(s))
    return (π/2) .* (s[:,2:end] * (1.0 ./ kk[2:end])) ./ (u_rms_sq) 
end

