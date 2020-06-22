using COVID
using Test

@testset "COVID.jl" begin
    params = ParameterSetting()
    res = run_model(params, .4, 20)
    @test typeof(res) == Array{Int64,2}
end
