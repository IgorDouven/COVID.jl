using COVID
using Test

@testset "COVID.jl" begin
    params = ParameterSetting()
    switch_pts = [(.3, 5), (.9, 5)]
    res1 = run_model(params, .4, 20)
    res2 = run_model(params, switch_pts)
    res3 = run_model(params, .03, 6, 20)
    @test typeof(res1) == Array{Int64,2}
    @test typeof(res2) == Array{Int64,2}
    @test typeof(res3) == Array{Int64,2}
end
