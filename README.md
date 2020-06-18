# COVID.jl

Julia code for agent-based model of COVID-19 outbreaks and possible intervention strategies.

## Installation
```julia
] add https://github.com/IgorDouven/COVID.jl
```
## Usage
To use default parameter setting:
```julia
using COVID

params = ParameterSetting()
```
This shows which fields can be altered. For instance, to change the maximum duration of infection, `max_tspan` (default = 24), to 20 days and the average household size, `hhs` (default = 4), to 6.5, do
```julia
params = ParameterSetting(max_tspan = 20, hhs = 6.5)
```
Similarly for the other parameters. NB: the output from calling `ParameterSetting()` also tells you the types of the various parameters; these have to be respected (so if, e.g., you want to set the average household size to 6, enter this as `ParameterSetting(hhs = 6.)`).

Then to let a model update for, say, 100 days, assuming that the probability that agents will stay home on any given day rather than visit one of their contacts equals .4, use
```julia
res = run_model(params, .4, 100)
```

Assuming this setting of this example, the following plots the infected and recovered (which includes the deceased) at all points in time:
```julia
cvd_plot(res)
```
![Results of updating the model for 100 time steps, showing the infected and recovered](./doc/IR.png)

If one would like to see a full SIRD output (so also plotting the susceptibles and showing separately the really recovered and the deceased), run
```julia
cvd_plot(res, sird=true)
```
![Same results, now showing also the susceptibles and separating the really recovered from the deceased](./doc/SIRD.png)

To study the effect of changing the probability of staying at home, one can do, for instance,
```julia
switch_pts = repeat([(.3, 10), (.9, 15), (.6, 10)], outer=3) 
res = run_model(params, switch_pts)
```
which updates the model during 10 steps assuming a probability of .3, followed by 15 steps assuming a probability of .9, followed by 10 steps assuming a probability of .6, and this repeated thrice. This yields the following result:
![Updating the model while switching the stay-at-home probability at various points in time](./doc/switch_SIRD.png)
