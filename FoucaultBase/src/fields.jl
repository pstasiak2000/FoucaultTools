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

export fieldkind
export dot!, cross!


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

struct FieldStyle <: Base.Broadcast.BroadcastStyle end

Base.BroadcastStyle(::Type{<:Field}) = FieldStyle()

Base.BroadcastStyle(
    ::FieldStyle,
    ::Base.Broadcast.DefaultArrayStyle
) = FieldStyle()
    
Base.BroadcastStyle(
    ::Base.Broadcast.DefaultArrayStyle,
    ::FieldStyle
) = FieldStyle()
    
Base.broadcastable(f::Field) = f

function Base.copy(bc::Base.Broadcast.Broadcasted{FieldStyle})
    data = Base.materialize(
        Base.Broadcast.Broadcasted(
            Base.Broadcast.DefaultArrayStyle{ndims(bc)}(),
            bc.f,
            bc.args
        )
    )

    Field(data)
end

Base.size(f::Field) = size(f.data)
Base.axes(f::Field) = axes(f.data)
Base.getindex(f::Field, I...) = getindex(f.data, I...)
Base.setindex!(f::Field, x, I...) = setindex!(f.data, x, I...)
Base.IndexStyle(::Type{<:Field{F, D, T, A}}) where {F,D, T, A} = IndexStyle(A)
Base.copy(f::Field) = Field(copy(f.data))


Base.zeros(::Type{<:FieldType{T}}, N::Vararg{Int,D}) where {T,D} = Field(zeros(T,N...))
Base.rand(::Type{<:FieldType{T}}, N::Vararg{Int,D}) where {T,D} = Field(rand(T,N...))
Base.similar(f::Field{F,D}) where {F,D} = Base.zeros(F,size(f)...)


fieldkind(::Field{F}) where {F} = F

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

Base.eltype(v::VectorField{F,C,D,T}) where {F,C,D,T} = T
Base.size(v::VectorField) = size(first(v.components))
Base.getindex(v::VectorField, i::Int) = v.components[i]
Base.setindex!(v::VectorField, x, I...) = setindex!(v.components, x, I...)
Base.length(v::VectorField{F,C}) where {F,C} = C
Base.ndims(v::VectorField) = length(v)
Base.copy(x::VectorField) = VectorField(map(copy, x.components))

Base.zeros(::Type{<:VectorField{F}},C::Int,M::Vararg{Int,D}) where {F,D} = VectorField(ntuple(_ -> zeros(F,M...), C))

Base.rand(::Type{<:VectorField{F}},C::Int,M::Vararg{Int,D}) where {F,D} = VectorField(ntuple(_ -> rand(F,M...), C))

Base.similar(v::VectorField{F,C}) where {F,C} = zeros(typeof(v),C,size(v)...)

fieldkind(::VectorField{F}) where {F} = F

check_size(x::VectorField,y::VectorField) = size(x) == size(y) || throw(FieldSizeMismatch(size(x),size(y)))

#Base.show(io::IO, ::MIME"text/plain", f::Field{F,D}) where {F,D} =
#    print(io, "Field{$F,$D} with size ", size(f))

"""
    Base.sum(v::VectorField)

Sum the components of the vector field piece-wise.
"""
Base.sum(v::VectorField) = sum(v.components)


####################################################################
#
#      Operations 
#
####################################################################

"""
    dot!(z::Field, x::VectorField, y::VectorField)

In-place variant of [`dot`](@ref) for vector fields.
"""
function dot!(z::Field{F,D}, x::VectorField{F,C,D}, y::VectorField{F,C,D}) where {F,C,D}
    @inbounds for i=1:C
        @. z += conj.(x[i]) * y[i]
    end
end

"""
    LinearAlgebra.dot(x::VectorField,y::VectorField)

Performs a piecewice dot product on two vector fields, and outputs a `Field` structure with the same spatial structure. 
"""
function LinearAlgebra.dot(x::VectorField{F,C,D},y::VectorField{F,C,D}) where {F,C,D}
    check_size(x,y)
    z = zeros(F,size(x)...)
    
    dot!(z,x,y)

    return z
end

"""
    FoucaultBase.cross!(z,x,y)

In-place variant of [`cross`](@ref) for vector fields.
"""
function cross!(z::VectorField{F,3,D},
                              x::VectorField{F,3,D}, 
                              y::VectorField{F,3,D}) where {F,D}
    @. z[1] = (x[2] * y[3]) - (x[3] * y[2])
    @. z[2] = (x[3] * y[1]) - (x[1] * y[3])
    @. z[3] = (x[1] * y[2]) - (x[2] * y[1])
    return nothing
end


"""
    LinearAlgebra.cross(x::VectorField, y::VectorField)

Performs a vector cross product on two vector fields, and outputs a `VectorField` structure type.
"""
function LinearAlgebra.cross(x::VectorField{F,3,D},y::VectorField{F,3,D}) where {F,D}
    check_size(x,y)
    z = zeros(typeof(x),3,size(x)...)
    cross!(z,x,y)
    return z
end

"""
    LinearAlgebra.norm(v::VectorField; p::Real=2)

Computes a p-norm of the vector field on the components.

The function defaults to the 2-norm.
"""
LinearAlgebra.norm(v::VectorField; p::Real=2) = sum(v .^ p)


