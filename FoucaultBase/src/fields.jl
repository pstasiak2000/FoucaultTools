### ====================================================== ##
### Field routines
### ====================================================== ##
#
#   This file consists of the definitions of fields and the routines involved in various computations directly from the 1D-3D.
# 
#

export Field

abstract type FieldType{T} end

struct RealField{T<:Real} <: FieldType{T} end
struct ComplexField{T<:Complex} <: FieldType{T} end


"""
    Field{D,A<:AbstractArray{T,D}} <: AbstractArray{T,D}

Main semantic wrapper used for fields.
"""
struct Field{F<:FieldType,D,T,A<:AbstractArray{T,D}} <: AbstractArray{T,D}
    data::A
end

Field(data::A) where {D,T<:Real,A<:AbstractArray{T,D}} = Field{RealField{T}, D, T, A}(data)
Field(data::A) where {D,T<:Complex,A<:AbstractArray{T,D}} = Field{ComplexField{T}, D, T, A}(data)

Base.size(f::Field) = size(f.data)
Base.getindex(f::Field, I...) = getindex(f.data, I...)
Base.setindex!(f::Field, x, I...) = setindex!(f.data, x, I...)
Base.IndexStyle(::Type{<:Field{F, D, T, A}}) where {F,D, T, A} = IndexStyle(A)

Base.zeros(::Type{<:FieldType{T}}, N::Vararg{Int,D}) where {T,D} = Field(zeros(T,N...))
Base.rand(::Type{<:FieldType{T}}, N::Vararg{Int,D}) where {T,D} = Field(rand(T,N...))


fieldtype(f::Field{F}) where {F} = F

component(f::Field{F,D}, i::Int) where {F,D} = selectdim(f, D, i) 




