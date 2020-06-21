using COVID
using Test

@testset "COVID.jl" begin
    pop, gr = create_model(100, .1, 4.)
    @test pop[end].pos == 100
    @test length(gr) == 1000
end
