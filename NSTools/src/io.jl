### ====================================================== ##
### Navier-Stokes I/O
### ====================================================== ##


function read_spectrum(path::AbstractString, file::Symbol)
    filename = get(HydroNS, file) do
        throw(ArgumentError("Unknown file key: $file"))
    end

    data = FoucaultBase.read_data(path, filename)
end

read_spectrum_time(path) = read_spectrum(path,:time)
