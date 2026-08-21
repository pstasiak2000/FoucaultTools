push!(LOAD_PATH,"../")

using FoucaultBase
using NSTools

using Test

path = "/home/piotr-stasiak/Codes/FoucaultTools/test/data"

tt = read_spectrum_time(path)

@testset verbose = true "NSTools Test set" begin

    @testset "HydroNS routines" begin
        @testset "Energy computation" begin
            E_kin = read_data(path,"energy.dat") 
            spU = read_spectrum(path,:energy)
            
            @test compute_energy(spU) ≈  E_kin[:,2]
        end

        @testset "Vorticity computation" begin
            vort = read_data(path,"vorticity.dat")
            spU = read_spectrum(path,:energy)
            
#            @test compute_enstrophy(spU) ≈ vort[:,2]
        end
    end

end
