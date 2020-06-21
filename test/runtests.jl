using COVID
using Test

@testset "COVID.jl" begin
    pop, gr = create_model(100, .1, 4.)
    pop[end].pos == 100
    length(gr) == 1000
end
