function cvd_plot(res::Array{Int64,2}; sird=false)
    l, _ = size(res)
    if sird == true
        plot(0:l - 1, res, 
            label = ["Susceptible" "Infected" "Recovered" "Deceased"], 
            xlabel = "Time",
            ylabel = "Count",
            linewidth = 5)
    else
        rec = res[:, 3] .+ res[:, 4]
        dat = hcat(res[:, 2], rec)    
        plot(0:l - 1, dat, 
            label = ["Infected" "Recovered"], 
            xlabel = "Time",
            ylabel = "Count",
            linewidth = 5)
    end
end

function cvd_plot(res::Array{Float64,2}; sird=false)
    l, _ = size(res)
    p1 = if sird == true
        plot(0:l - 1, res[:, 1:4], 
            label = ["Susceptible" "Infected" "Recovered" "Deceased"],
            ylabel = "Count",
            linewidth = 5)
    else
        rec = res[:, 3] .+ res[:, 4]
        dat = hcat(res[:, 2], rec)    
        plot(0:l - 1, dat, 
            label = ["Infected" "Recovered"], 
            xlabel = "Time",
            ylabel = "Count",
            linewidth = 5)
    end
    p2 = plot(0:l - 1, res[:, 5], legend = false, xlabel = "Time", ylabel = "Weight", linewidth = 5)
    return plot(p1, p2, layout = (2, 1))
end
