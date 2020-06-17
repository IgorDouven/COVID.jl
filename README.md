# COVID.jl

Julia code for agent-based model of COVID-19 outbreaks and possible intervention strategies.

## Installation
```julia
] add https://github.com/IgorDouven/COVID.jl
```
## Usage
To use default model parameters:
```julia
mod = Model()
```
This shows which fields can be altered. For instance, to change the maximum duration of infection, `max_tspan`, to 20 days, do
```julia
mod = Model(max_tspan = 20)
```
Similarly for the other parameters.

Then to let a model update for, say, 100 days, assuming that the probability that agents will stay home on any given day rather than visit one of their contacts equals .4, use
```julia
res = sim_run(mod, .4, 100)
```

The following plots the infected and recovered (including deceased) at any point in time:
```julia
cvd_plot(res)
```
If one would like to see a full SIRD output (so also plotting the susceptibles and showing separately the really recovered and the deceased), run
```julia
cvd_plot(res, sird=true)
```
