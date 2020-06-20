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

# the following is Deb's fast non dominating sorting algorithm
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
