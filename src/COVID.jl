"""
    COVID

An agent-based approach to modeling COVID-19 outbreaks and possible intervention strategies
"""
module COVID

export
    ParameterSetting,
    run_model,
    run_evo,
    cvd_plot

using Distributed
using LinearAlgebra
using Distributions
using Distances
using DataFrames
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
