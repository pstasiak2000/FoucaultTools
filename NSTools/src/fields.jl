### ====================================================== ##
### Field routines
### ====================================================== ##
#
#   This file consists of the definitions of fields and the routines involved in various computations directly from the 1D-3D velocity and other fields derived from foucault.
# 
#

abstract type FieldType end

struct RealField <: FieldType end
struct ComplexField <: FieldType end


"""
    Field{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}

Main semantic wrapper used for fields.
"""
struct Field{F<:FieldType,T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
    data::A
end

Field(data::A) where {T<:Real,N,A<:AbstractArray{T,N}} = Field{RealField, T, N, A}(data)
Field(data::A) where {T<:Complex,N,A<:AbstractArray{T,N}} = Field{ComplexField, T, N, A}(data)

Base.size(f::Field) = size(f.data)
Base.getindex(f::Field, I...) = getindex(f.data, I...)
Base.setindex!(f::Field, x, I...) = setindex!(f.data, x, I...)
Base.IndexStyle(::Type{<:Field{F,T,N,A}}) where {F,T,N,A} = IndexStyle(A)

fieldtype(f::Field{F,T,N,A}) where {F,T,N,A} = F
