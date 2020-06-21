struct Proposal
    W::Int
    S::Float64
end

rand_proposal() = Proposal(rand(5:15), rand(Uniform(0, .1)))

function crowdingDistance(arr::Array{Float64,2})
    rs, cs = size(arr)
    c_sum = zeros(rs)
    for k ∈ 1:cs
        a = arr[:, k]
        cd = Vector{Float64}(undef, rs)
        cd[1] = cd[end] = Inf
        sc = sort(a)
        ind = [ sort(collect(enumerate(a)), by = x -> x[2])[i][1] for i in 1:rs ]
        dist = euclidean(extrema(a)...)
        @inbounds for i in 2:rs - 1
            cd[i] = (sc[i + 1] - sc[i - 1])/dist
        end
        c_sum += cd[ind, :]
    end
    return c_sum
end

function dominates(x, y)
    strict_inequality_found = false
    for i ∈ eachindex(x)
        y[i] < x[i] && return false
        strict_inequality_found |= x[i] < y[i]
    end
    return strict_inequality_found
end

function fast_nds(arr)
    sz = size(arr, 2)
    fronts = Vector{Int64}[]
    ind = collect(axes(arr, 1))
    a = SVector{sz}.(eachrow(arr))
    N = length(a)
    ii = 1
    while ii <= N
        start = ii
        @inbounds for i in ii:N
            for j in start:N
                dominates(a[j], a[i]) && @goto l1
            end
            ind[ii], ind[i] = ind[i], ind[ii]
            a[ii], a[i] = a[i], a[ii]
            ii += 1
            @label l1
        end
        push!(fronts, ind[start:ii - 1])
    end
    return fronts
end

function mix(p::Proposal, q::Proposal)
    w = p.W < q.W ? rand(p.W:q.W) : rand(q.W:p.W)
    s = p.S < q.S ? rand(Uniform(p.S, q.S)) : rand(Uniform(q.S, p.S))
    return Proposal(w, s)
end

function crossover(v::Array{Proposal,1})
    parent_pop_size = length(v)
    child_arr = Array{Proposal,1}(undef, parent_pop_size)
    @inbounds for i in 1:parent_pop_size
        x, y = sample(1:parent_pop_size, 2, replace = false)
	child_arr[i] = mix(v[x], v[y])
    end
    return child_arr
end

function new_gen(v::Array{Proposal,1}; prob = 0.05)
    pop = vcat(v, crossover(v))
    s = sample(1:length(pop), floor(Int, 100*prob), replace = false)
    @inbounds for i ∈ s
	pop[i] = rand_proposal()
    end
    return pop
end

function intervene(ps::ParameterSetting, p::Proposal, numb_updates::Int; init_w = .5)
    @unpack max_tspan, pr_quick_rec, trans_rate, drop, pr_death, N, λ, hhs = ps
    ppl, grph = create_model(N, λ, hhs)
    window, scale = p.W, p.S
    daily_inf = Int[]
    weights = Float64[]
    cfs = Float64[]
    
    @inbounds for i in 1:window
        model_update!(ppl, grph, init_w, max_tspan, pr_quick_rec, trans_rate, drop, pr_death)
        push!(daily_inf, sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :I)))
        push!(weights, init_w)
    end

    pq = .0
    @inbounds for i in 1:numb_updates - window
        x = [ones(window) 1:window]
        y = daily_inf[end - (window - 1):end]
        cf = x \ y
        pq = pq + (cf[2] * scale)
        w = rich_scaled(pq)
        model_update!(ppl, grph, w, max_tspan, pr_quick_rec, trans_rate, drop, pr_death)
        push!(daily_inf, sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :I)))
        push!(weights, w)
    end

    mw = mean(weights)
    undercap = sum(max.((daily_inf ./ 100) .- (length(ppl) / 1000), 0))
    tot_dec = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :D))

    return tot_dec, mw, undercap
end

function sim_fnc_nsga(ps::ParameterSetting, p::Proposal, numb_updates::Int; n_test=5)
    dec = 0
    uc = 0.
    wgts = 0.

    i = 1
    while i <= n_test
        res = intervene(ps, p, numb_updates)
        if first(res) < 4 # epidemic didn't get off the ground
            i -= 1
        else
            dec += first(res)
            uc += last(res)
            wgts += res[2]
        end
        i += 1
    end
    return dec / n_test, uc / n_test, wgts / n_test
end

function select_best(ps::ParameterSetting, arr::Array{Proposal,1}, numb_updates::Int, parent_pop_size::Int)
    scores = Matrix(DataFrame(pmap(i -> sim_fnc_nsga(ps, arr[i], numb_updates), 1:length(arr))))
    cdist = crowdingDistance(scores)
    ranking = fast_nds(scores)
    rnkInd = Vector{Int64}(undef, 2*parent_pop_size)
    @inbounds for i in 1:length(ranking)
        rnkInd[ranking[i]] .= i
    end
    ns = hcat(rnkInd, cdist[:], 1:2*parent_pop_size)
    p2 = sortperm(view(ns, :, 2))
    s1 = ns[p2, :]
    p1 = sortperm(view(s1, :, 1))
    out = Int.(s1[p1, 3][1:parent_pop_size])
    return arr[out], scores[out, :]
end

extract_pams(p::Proposal) = (p.W, p.S)

function run_evo(ps::ParameterSetting, numb_updates::Int, numb_gen::Int, parent_pop_size::Int; full=false)
    pms_ar = Array{Tuple{Int64,Float64},1}[]
    start_pop = [ rand_proposal() for _ in 1:parent_pop_size ]
    push!(pms_ar, extract_pams.(start_pop))
    scrs = Array{Float64,3}(undef, parent_pop_size, 3, numb_gen)
    for i in 1:numb_gen
        new_pop = [ Proposal(pms_ar[i][j]...) for j in 1:parent_pop_size ]
        out_pop, out_scrs = select_best(ps, new_gen(new_pop), numb_updates, parent_pop_size)
        push!(pms_ar, extract_pams.(out_pop))
        scrs[:, :, i] = out_scrs
    end
	
    pf = fast_nds(scrs[:, :, end])
    pareto_front = pms_ar[end][pf[1]]

    return full == true ? pms_ar, scrs, pareto_front : pareto_front
end
