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
    VectorField{F<:FieldType, C, D, T, A}

Wrapper used for vector fields.
"""
struct VectorField{F,C,D,T,A<:AbstractArray{T,D}}
    components::NTuple{C,Field{F,D,T,A}} 
end

const RealVectorField{T<:Real} = VectorField{RealField{T}}
const ComplexVectorField{T<:Complex} = VectorField{ComplexField{T}} 

const RealVectorField64 = VectorField{RealField64}
const ComplexVectorField64 = VectorField{ComplexField64}

VectorField(components::NTuple{C,Field{F,D,T,A}}) where {C,D,F,T,A} = VectorField{F,C,D,T,A}(components)
VectorField(components::NTuple{C,AbstractArray{T,D}}) where {C,D,T} = VectorField(Field.(components))

VectorField(components::Vararg{Field}) = VectorField(components)






struct VectorFieldStyle <: Base.Broadcast.BroadcastStyle end

Base.BroadcastStyle(::Type{<:VectorField}) = VectorFieldStyle()
Base.broadcastable(x::VectorField) = x

Base.BroadcastStyle(
    ::VectorFieldStyle,
    ::Base.Broadcast.DefaultArrayStyle
) = VectorFieldStyle()

Base.BroadcastStyle(
    ::Base.Broadcast.DefaultArrayStyle,
    ::VectorFieldStyle
) = VectorFieldStyle()

_broadcast_component(x::VectorField, i) = x[i]
_broadcast_component(x, i) = x

function _broadcast_component(
    bc::Base.Broadcast.Broadcasted,
    i
)
    Base.Broadcast.broadcasted(
        bc.f,
        (_broadcast_component(arg, i) for arg in bc.args)...
    )
end

function Base.copy(
    bc::Base.Broadcast.Broadcasted{VectorFieldStyle}
)
    vf = _find_vectorfield(bc)

    C = length(vf.components)

    components = ntuple(Val(C)) do i
        Base.materialize(_broadcast_component(bc, i))
    end

    VectorField(components)
end

function _find_vectorfield(bc::Base.Broadcast.Broadcasted)
    for arg in bc.args
        if arg isa VectorField
            return arg
        elseif arg isa Base.Broadcast.Broadcasted
            return _find_vectorfield(arg)
        end
    end

    error("No VectorField found in broadcast expression")
end

_find_vectorfield(x::VectorField) = x





Base.size(v::VectorField) = size(first(v.components))
Base.getindex(v::VectorField, i::Int) = v.components[i]
Base.setindex!(v::VectorField, x, I...) = setindex!(v.components, x, I...)
Base.length(v::VectorField{F,C}) where {F,C} = C

#
#function Base.:+(x::VectorField{F,C,D},
#        y::VectorField{F,C,D}
#       ) where {F,C,D}
#    VectorField(ntuple(i -> x[i] + y[i],Val(C)))  
#end
#
#function Base.:*(c::Number, x::VectorField{F,C,D}) where {F,C,D}
#    VectorField(
#        ntuple(i -> c .* x.components[i], Val(C))
#    )
#end
#
#Base.:*(x::VectorField, c::Number) = c * x
#
#function Base.broadcast(
#        ::typeof(*),
#        x::VectorField{F,C,D},
#        y::VectorField{F,C,D}
#    ) where {F,C,D}
#    
#    VectorField(ntuple(i -> x[i] .* y[i],Val(C)))
#end

Base.zeros(::Type{<:VectorField{F}},C::Int,M::Vararg{Int,D}) where {F,D} = VectorField(ntuple(_ -> zeros(F,M...), C))
Base.rand(::Type{<:VectorField{F}},C::Int,M::Vararg{Int,D}) where {F,D} = VectorField(ntuple(_ -> rand(F,M...), C))

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



