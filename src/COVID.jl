"""
    COVID

An agent-based approach to modeling COVID-19 outbreaks and possible intervention strategies
"""
module COVID

export
	Model,
    create_model,
    model_update!,
    abm_run,
    sim_run,
    cvd_plot

using Distributions
using Random
using StatsBase
using LightGraphs
using Parameters
using Plots

include("COVID_FNCS.jl")
include("COVID_PLOT.jl")

end