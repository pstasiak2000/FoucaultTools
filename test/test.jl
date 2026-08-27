push!(LOAD_PATH,"../")

using FoucaultBase
using NSTools

using Test

#path = "/home/piotr-stasiak/Codes/FoucaultTools/test/data"
path = joinpath(@__DIR__,"data")

tt = read_spectrum_time(path)

@testset verbose = true "FoucaultBase Test set" begin
    
    @testset verbose = true "Testing Fields" begin
        N = 64
        x = rand(RealField64,N,N,N)
        y = rand(ComplexField64,N,N,N)

        @testset verbose = true "Field Type stability" begin
            @test typeof(x + x) <: typeof(x)
            @test typeof(2.0 .* x) <: typeof(x)
            @test typeof(x .* x) <: typeof(x)

            @test typeof(y + y) <: typeof(y)
            @test typeof(2.0 .* y) <: typeof(y)
            @test typeof(y .* conj.(y)) <: typeof(y)
            
            @test typeof(real.(y)) <: Field{RealField64}
        end
    end

    @testset verbose = true "Testing VectorFields" begin
        N = 64
        v = rand(RealVectorField64,3,N,N,N)
        w = rand(ComplexVectorField64,3,N,N,N)

        a = similar(v)
        b = similar(w)
        
        f = rand(fieldkind(v),size(v)...)
        g = rand(fieldkind(w),size(w)...)


        @testset "VectorField Type Stability" begin
            @test typeof(v .+ v) <: typeof(v)
            @test typeof(2.0 .* v) <: typeof(v)
            @test typeof(v .* v) <: typeof(v)

            @test typeof(w .+ w) <: typeof(w)
            @test typeof(2.0 .* w) <: typeof(w)
            @test typeof(w .* conj.(w)) <: typeof(w)

            @test typeof(real.(w)) <: VectorField{RealField64}
        end
           

        @testset "Vector field op. allocs." begin
            @test (@allocated dot!(f,v,v)) == 0
            @test (@allocated dot!(g,w,w)) == 0
            
            @test (@allocated cross!(a,v,v)) == 0
            @test (@allocated cross!(b,w,w)) == 0
        end
    end

end


#@testset verbose = true "NSTools Test set" begin
#
#    @testset "HydroNS routines" begin
#        @testset "Energy computation" begin
#            E_kin = read_data(path,"energy.dat") 
#            spU = read_spectrum(path,:energy)
#            
#            @test compute_energy(spU) ≈  E_kin[:,2]
#        end
#
#        @testset "Vorticity computation" begin
#            vort = read_data(path,"vorticity.dat")
#            spU = read_spectrum(path,:energy)
#            
##            @test compute_enstrophy(spU) ≈ vort[:,2]
#        end
#    end
#
#end
