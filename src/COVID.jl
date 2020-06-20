"""
    COVID

An agent-based approach to modeling COVID-19 outbreaks and possible intervention strategies
"""
module COVID

export
    ParameterSetting,
    create_model,
    model_update!,
    abm_run,
    run_model,
    run_flex,
    cvd_plot

using Distributions
using Random
using StatsBase
using StaticArrays
using LightGraphs
using Parameters
using Plots
using Pipe

include("COVID_FNCS.jl")
include("COVID_FLEX.jl")
include("nsga2.jl")
include("COVID_PLOT.jl")

end
