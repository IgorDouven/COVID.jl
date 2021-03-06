rich_scaled(x::Float64) = 0.95 / (1 + exp(-x))

function run_model(ps::ParameterSetting, scale::Float64, window::Int, numb_updates::Int; init_w=0.5)
    @unpack max_tspan, pr_quick_rec, trans_rate, drop, pr_death, N, λ, hhs = ps
    ppl, grph = create_model(N, λ, hhs)
    s = vcat(dropdims(sum(retrieve_data(ppl), dims=1), dims=1), init_w)
    daily_sus = Vector{Int64}(undef, numb_updates)
    daily_inf = Vector{Int64}(undef, numb_updates)
    daily_rec = Vector{Int64}(undef, numb_updates)
    daily_dead = Vector{Int64}(undef, numb_updates)
    wgts = Vector{Float64}(undef, numb_updates)
    
    @inbounds for i in 1:window
        model_update!(ppl, grph, init_w, max_tspan, pr_quick_rec, trans_rate, drop, pr_death)
        daily_sus[i] = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :S))
        daily_inf[i] = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :I))
        daily_rec[i] = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :R))
        daily_dead[i] = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :D))
        wgts[i] = init_w
    end

    ps = 0.
    @inbounds for i in window + 1:numb_updates
        x = [ones(window) 1:window]
        y = daily_inf[i - window:i - 1]
        _, cf = x \ y
        ps = ps + (cf * scale)
        w = rich_scaled(ps)
        model_update!(ppl, grph, w, max_tspan, pr_quick_rec, trans_rate, drop, pr_death)
        daily_sus[i] = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :S))
        daily_inf[i] = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :I))
        daily_rec[i] = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :R))
        daily_dead[i] = sum(Int.([ ppl[j].status for j in 1:length(ppl) ] .== :D))
        wgts[i] = w
    end
    return vcat(s', hcat(daily_sus, daily_inf, daily_rec, daily_dead, wgts))
end
