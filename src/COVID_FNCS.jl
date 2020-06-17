# default parameter settings justified by what has been reported about coronavirus disease 2019 (COVID-19) up till beginning of June, 2020
@with_kw struct Model
    max_tspan::Int = 24 # maximum duration of infection
    pr_quick_rec::Float64 = .5 # half the infected recover in `max_tspan`/2 days
    trans_rate::Float64 = .04 # transmission rate
    drop::Int = 4 # by how much the infectiousness drops if the infection last longer than half of `max_tspan`
    pr_death::Float64 = .1 # probability of death
    N::Int = 1000 # number of households
    λ::Float64 = .0125 # probability that two nodes are connected
end

mutable struct Agent
    pos::Int
    tspan::Int
    status::Symbol
end

function create_model(N::Int, λ::Float64)
	g = LightGraphs.SimpleGraph(LightGraphs.Generators.ErdosRenyi(N, round(Int, binomial(N, 2)*λ)))	
	pop = Agent[]
	for i in 1:N
		r = 1 + rand(Poisson(3)) 
    		for j in 1:r
			push!(pop, Agent(i, 0, :S))
		end
	end
	links = [ length(neighbors(g, i)) for i in 1:N ]
	quasi_isolated = collect(1:N)[links .< 2]
    inds = [ pop[i].pos ∉ quasi_isolated for i in 1:length(pop) ]
    ag_zeros = sample(collect(1:length(pop))[inds], length(pop) ÷ 100, replace=false) # around 1 percent of population is infected
    for i in ag_zeros
        pop[i] = Agent(pop[i].pos, 1, :I)
    end
	return pop, g
end

function meet(p::Array{Agent,1}, 
              g::LightGraphs.SimpleGraphsCore.SimpleGraph{Int64}, 
              a::Agent, 
              w::Float64,
              max_tspan::Int,
              pr_quick_rec::Float64)
    ts = max_tspan * pr_quick_rec
    if a.status == :S || (0 < a.tspan <= ts)
        nbs = neighbors(g, a.pos) # list all neighbors
        weight_lst = Int[]
        @inbounds for i in nbs
            p_sel_ind = [ p[j].pos == i for j in 1:length(p) ]
            p_sel = p[p_sel_ind]
            sel_ag = [ p_sel[j].status == :S || (0 < p_sel[j].tspan <= ts) for j in 1:length(p_sel) ]
            push!(weight_lst, length(sel_ag))
        end
        s = ifelse(isempty(weight_lst), a.pos, sample(nbs, Weights(weight_lst)))
        v = ifelse(rand(Bernoulli(w)), a.pos, s)
        return v
    end
end

function contagion!(p::Array{Agent,1}, 
                    g::LightGraphs.SimpleGraphsCore.SimpleGraph{Int64}, 
                    w::Float64,
                    max_tspan::Int,
                    pr_quick_rec::Float64,
                    trans_rate::Float64,
                    drop::Int)
    ts = max_tspan * pr_quick_rec
    s = shuffle(1:length(p))
    @inbounds for a in s
        m = meet(p, g, p[a], w, max_tspan, pr_quick_rec)
        hh_sel_ind = [ p[i].pos == m for i in 1:length(p) ]
        hh_sel = p[hh_sel_ind]
        group = vcat(hh_sel, p[a])
        numb_inf_short = sum([ group[i].status == :I for i in 1:length(group) ] ∩ 
            [ group[i].tspan <= ts for i in 1:length(group) ])
        numb_inf_long = sum([ group[i].status == :I for i in 1:length(group) ] ∩ 
            [ group[i].tspan > ts for i in 1:length(group) ])
        while (numb_inf_short = numb_inf_short - 1) >= 0
            @inbounds for i in 1:length(group)
                if group[i].status == :S && rand(Bernoulli(trans_rate))
                    group[i].status = :I
                end
            end
        end
        while (numb_inf_long = numb_inf_long - 1) >= 0
            @inbounds for i in 1:length(group)
                if group[i].status == :S && rand(Bernoulli(trans_rate/drop))
                    group[i].status = :I
                end
            end
        end
    end
end

function update_infection_duration!(p::Array{Agent,1})
    ind = [ p[i].status == :I for i in 1:length(p) ]
    p = p[ind]
    @inbounds for i in 1:length(p)
        p[i].tspan += 1
    end
end

function recover_or_not!(p::Array{Agent,1},
                         max_tspan::Int,
                         pr_quick_rec::Float64,
                         pr_death::Float64)
    ind_short = [ p[i].tspan == max_tspan * pr_quick_rec for i in 1:length(p) ]
    ind_long = [ p[i].tspan == max_tspan for i in 1:length(p) ]
    p_short = p[ind_short]
    p_long = p[ind_long]
    @inbounds for i in 1:length(p_short)
        if rand(Bernoulli(pr_quick_rec))
            p_short[i].status = :R
            p_short[i].tspan = 0
        end
    end
    @inbounds for i in 1:length(p_long)
        if rand(Bernoulli(pr_death))
            p_long[i].status = :D
        else
            p_long[i].status = :R
        end
        p_long[i].tspan = 0
    end
end

function retrieve_data(p::Array{Agent,1})
	res = Array{Int64,2}(undef, p[end].pos, 4)
	@inbounds for n in 1:p[end].pos
        p_sel_ind = [ p[i].pos == n for i in 1:length(p) ]
        p_sel = p[p_sel_ind]
        household_stat = [ p_sel[i].status for i in 1:length(p_sel) ]
        res[n, :] = [ count(i->i==X, household_stat) for X in [:S, :I, :R, :D] ]
    end
    return res
end

function model_update!(p::Array{Agent,1}, 
                       g::LightGraphs.SimpleGraphsCore.SimpleGraph{Int64}, 
                       w::Float64, 
                       max_tspan::Int, 
                       pr_quick_rec::Float64, 
                       trans_rate::Float64,
                       drop::Int,
                       pr_death::Float64)
    contagion!(p, g, w, max_tspan, pr_quick_rec, trans_rate, drop)
    update_infection_duration!(p)
    recover_or_not!(p, max_tspan, pr_quick_rec, pr_death)
    return retrieve_data(p)
end

function abm_run(p::Array{Agent,1}, 
                 g::LightGraphs.SimpleGraphsCore.SimpleGraph{Int64}, 
                 w::Float64, 
                 mt::Int,
                 pqr::Float64,
                 tr::Float64,
                 d::Int,
                 prd::Float64,
                 numb_updates::Int)
    p_dc = deepcopy(p)
    res_array = [ model_update!(p_dc, g, w, mt, pqr, tr, d, prd) for _ in 1:numb_updates ]
    res_totals = [ sum(res_array[i], dims=1) for i in 1:numb_updates ]
    out = reduce(vcat, res_totals)
    return out
end

function sim_run(m::Model, w::Float64, numb_updates::Int)
    @unpack max_tspan, pr_quick_rec, trans_rate, drop, pr_death, N, λ = m
    population, graph = create_model(N, λ)
    out = abm_run(population, graph, w, max_tspan, pr_quick_rec, trans_rate, drop, pr_death, numb_updates)
    s = sum(retrieve_data(population), dims=1)
    return vcat(s, out)
end
