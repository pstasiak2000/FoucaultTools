### ====================================================== ##
### Field routines
### ====================================================== ##
#
#   This file consists of the definitions of fields and the routines involved in various computations directly from the 1D-3D.
# 
#

export Field
export VectorField

export RealField, ComplexField
export RealField64, ComplexField64

export RealVectorField, ComplexVectorField 
export RealVectorField64, ComplexVectorField64

export dot!

abstract type FieldType{T} end

struct RealField{T<:Real} <: FieldType{T} end
struct ComplexField{T<:Complex} <: FieldType{T} end

const RealField64 = RealField{Float64}
const ComplexField64 = ComplexField{ComplexF64}

struct FieldSizeMismatch <: Exception
    xsize
    ysize
end

function Base.showerror(io::IO, e::FieldSizeMismatch) 
    println("Field size mismatch between vector field components: Both fields must have the same size")
    println("Field 1: $(e.xsize)")
    println("Field 2: $(e.ysize)")
end


####################################################################
#
#       (Scalar) Field definitions
#
####################################################################

"""
    Field{D,A<:AbstractArray{T,D}} <: AbstractArray{T,D}

Main semantic wrapper used for a scalar fields.
"""
struct Field{F<:FieldType,D,T,A<:AbstractArray{T,D}} <: AbstractArray{T,D}
    data::A
end

Field(data::A) where {D,T<:Real,A<:AbstractArray{T,D}} = Field{RealField{T}, D, T, A}(data)
Field(data::A) where {D,T<:Complex,A<:AbstractArray{T,D}} = Field{ComplexField{T}, D, T, A}(data)

Base.size(f::Field) = size(f.data)
Base.axes(f::Field) = axes(f.data)
Base.getindex(f::Field, I...) = getindex(f.data, I...)
Base.setindex!(f::Field, x, I...) = setindex!(f.data, x, I...)
Base.IndexStyle(::Type{<:Field{F, D, T, A}}) where {F,D, T, A} = IndexStyle(A)

Base.zeros(::Type{<:FieldType{T}}, N::Vararg{Int,D}) where {T,D} = Field(zeros(T,N...))
Base.rand(::Type{<:FieldType{T}}, N::Vararg{Int,D}) where {T,D} = Field(rand(T,N...))

fieldtype(f::Field{F}) where {F} = F

####################################################################
#
#       Vector Field definitions
#
####################################################################
"""
    VectorField{F<:FieldType, D, T}

Wrapper used for vector fields.
"""
struct VectorField{F,D,T,A<:AbstractArray{T,D}}
    components::NTuple{D,Field{F,D,T,A}} 
end

const RealVectorField{T<:Real} = VectorField{RealField{T}}
const ComplexVectorField{T<:Complex} = VectorField{ComplexField{T}} 

const RealVectorField64 = VectorField{RealField64}
const ComplexVectorField64 = VectorField{ComplexField64}

VectorField(components::NTuple{D,Field{F,D,T,A}}) where {D,F,T,A} = VectorField{F,D,T,A}(components)
VectorField(components::NTuple{D,AbstractArray{T,D}}) where {D,T} = VectorField(Field.(components))

VectorField(components::Vararg{Field}) = VectorField(components)

Base.size(v::VectorField) = size(first(v.components))
Base.getindex(v::VectorField, i::Int) = v.components[i]
Base.setindex!(v::VectorField, x, I...) = setindex!(v.components, x, I...)

Base.zeros(::Type{<:VectorField{F}},M::Vararg{Int,D}) where {F,D} = VectorField(ntuple(_ -> zeros(F,M...),Val(D)))
Base.rand(::Type{<:VectorField{F}},M::Vararg{Int,D}) where {F,D} = VectorField(ntuple(_ -> rand(F,M...),Val(D)))

fieldtype(v::VectorField{F}) where {F} = F

check_size(x::VectorField,y::VectorField) = size(x) == size(y) || throw(FieldSizeMismatch(size(x),size(y)))


####################################################################
#
#      Operations 
#
####################################################################

"""
    dot!(z::Field, x::VectorField, y::VectorField)

In-place variant of [`dot`](@ref) for vector fields.
"""
function dot!(z::Field{F,D}, x::VectorField{F,D}, y::VectorField{F,D}) where {F,D}
    @inbounds for i=1:D
        @. z += conj.(x[i]) * y[i]
    end
end

"""
    LinearAlgebra.dot(x::VectorField,y::VectorField})

Performs a piecewice dot product on two vector fields, and outputs a `Field` structure with the same spatial structure. 
"""
function LinearAlgebra.dot(x::VectorField{F,D},y::VectorField{F,D}) where {F,D}
    check_size(x,y)
    z = zeros(F,size(x)...)
    
    dot!(z,x,y)

    return z
end



