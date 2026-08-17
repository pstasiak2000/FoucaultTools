module FoucaultBase

using DelimitedFiles

"""
    read_data(dir, file, ::Type{T})

Read a formatted data file `file` from `dir` into a matrix of data type T.
"""
function read_data(
                    dir::AbstractString,
                    file::AbstractString,
                    ::Type{T}=Float64,
                   ) where {T}
    
    raw = readdlm(joinpath(dir,file))
    
    ### Handles reading of Fortran exponent notation and converts
    ### all string data into numeric type T (by default Float64)
    data = map(raw) do x
        if x isa AbstractString
            parse(T, replace(x, "D" => "E"))
        else
            T(x)
        end
    end
    return data
end

end # module FoucaultBase
