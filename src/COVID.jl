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
    cvd_plot

using Distributions
using Random
using StatsBase
using LightGraphs
using Parameters
using Plots
using Pipe

include("COVID_FNCS.jl")
include("COVID_PLOT.jl")

end
