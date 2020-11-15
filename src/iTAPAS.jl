using DataFrames
using CSV
using Random
using Calculus
using Printf
using StatsBase
using Dates
Random.seed!(1403)

@doc "
    traffic_assignment(;networkName, tol=1e-5, maxIters=20, maxRunTime=600, log=:on)

    improved Traffic Assignment by Paired Alternative Segments (iTAPAS) algorithm
    for static multi-class traffic assignment problem with generalized link cost
    function.

    ## Generalized link cost function: c(v) = Σ ℿₚ * p(v) * t
    -   c(v): generalized link cost for link ij
    -   t   : travel time on link ij
    -   v   : travel speed on link ij
    -   ℿₚ  : cost of parameter p
    -   p(v): parameter p as a polynomial function of v

    ## Required properties of the generalized cost function
    -   Strictly positive
    -   Monotonically non-decreasing
    -   COntinuously differentiable

    ## Arguments
    -   networkName : network from the repository https://github.com/anmol1104/TrafficAssignment
    -   tol         : tolerance level for relative gap convergence
    -   maxIters    : maximum number of iterations
    -   maxRunTime  : maximum wall clock run time (s)
    -   log         : presents results for every iteration if log is on

    ## DataFiles (available at: https://github.com/anmol1104/TrafficAssignment)
    -   cost    : Enlists cost (ℿₚ) for all the parameters (p) of the generalized cost function
    -   coef    : Enlists coefficients of p(v) for all the parameters (p) of the generalized cost function
    -   class   : Enlists the relevant subset of parameters for the generalized cost function for each class
    -   network : Details the topology of the network
    -   demand  : Enlists OD pairs and corresponding demand

    ## IO Units
    -   length  : miles
    -   time    : hour
    -   volume  : litre
    -   mass    : kg
    -   cost    : \$
"
function traffic_assignment(;networkName, tol=1e-5, maxIters=20, maxRunTime=600, log=:on)
    println()
    printstyled("\niTAPAS Algorithm", color=:blue)

    N   = Int64[]                               # Nodes
    A   = Array{Int64,1}[]                      # Arcs as adjacency list
    M   = Int64[]                               # Modes/classes
    Vᵢⱼ = Array{Float64,1}[]                    # Link volume capcity
    dᵢⱼ = Array{Float64,1}[]                    # Link length
    tᵢⱼ = Array{Float64,1}[]                    # Link free flow travel time
    αᵢⱼ = Array{Float64,1}[]                    # BPR parameters
    βᵢⱼ = Array{Float64,1}[]                    # BPR parameters
    R   = Int64[]                               # Origins
    Mᵣ  = Dict{Int64,Int64}()                   # Mode/class assosciated with every origin
    Sᵣ  = Dict{Int64,Array{Int64,1}}()          # Destinations for every origin
    qᵣ  = Dict{Tuple{Int64,Int64},Float64}()    # Demand between OD pairs
    η   = Array{Float64,1}[]                    # Effective coeffficients for the polynomial c(v) = Σ ℿₚ * p(v) * t
    ϕ   = Array{Int64,1}[]


    function build()
        # cost file
        costFile = "src\\Network\\$networkName\\cost.csv"
        csv₁ = CSV.File(costFile, types=[String, Float64])
        df₁ = DataFrame(csv₁)
        parameters = [df₁[i,1] for i in 1:nrow(df₁)]::Array{String,1}
        ℿ = df₁[!,2]::Array{Float64,1}

        # coef file
        coefFile = "src\\Network\\$networkName\\coef.csv"
        csv₂ = CSV.File(coefFile)
        df₂ = DataFrame(csv₂)
        γ = [[df₂[i,j] for j in 2:ncol(df₂)] for i in 1:length(parameters)]::Array{Array{Float64,1},1}

        # criteria file
        clssFile = "src\\Network\\$networkName\\class.csv"
        csv₃ = CSV.File(clssFile, types=[Int64, String])
        df₃ = DataFrame(csv₃)
        criteria = [split(df₃[!,2][i], ", ") for i in 1:nrow(df₃)]::Array{Array{SubString{String},1},1}
        for m in df₃[!,1]::Array{Int64,1}
            append!(M,m)
            push!(η, [])
            for i in 1:(ncol(df₂)-1)
                append!(η[end], 0.0)
                for j in 1:length(parameters)
                    p = parameters[j]
                    if p ∈ criteria[m] η[m][i] += γ[j][i] * ℿ[j] end
                end
            end
        end

        # network file
        ntwkFile = "src\\Network\\$networkName\\network.csv"
        csv₄ = CSV.File(ntwkFile, types=[Int64, Int64, Float64, Float64, Float64, Float64, Float64])
        df₄ = DataFrame(csv₄)
        head = df₄[!, 1]::Array{Int64,1}
        tail = df₄[!, 2]::Array{Int64,1}
        linkcapacity = df₄[!, 3]::Array{Float64,1}
        linklength = df₄[!, 4]::Array{Float64,1}
        linkfft = df₄[!, 5]::Array{Float64,1}
        alpha = df₄[!, 6]::Array{Float64,1}
        beta = df₄[!, 7]::Array{Float64,1}
        n = max(maximum(head), maximum(tail))
        for i in 1:n
            append!(N, i)
            push!(A, [])
            push!(Vᵢⱼ, [])
            push!(dᵢⱼ, [])
            push!(tᵢⱼ, [])
            push!(αᵢⱼ, [])
            push!(βᵢⱼ, [])
            push!(ϕ, [])
        end
        for i in 1:length(head)
            append!(A[head[i]], tail[i])
            append!(Vᵢⱼ[head[i]], linkcapacity[i])
            append!(dᵢⱼ[head[i]], linklength[i])
            append!(tᵢⱼ[head[i]], linkfft[i])
            append!(αᵢⱼ[head[i]], alpha[i])
            append!(βᵢⱼ[head[i]], beta[i])
            append!(ϕ[head[i]], 1)
        end

        # demand file
        dmndFile = "src\\Network\\$networkName\\demand.csv"
        csv₅ = CSV.File(dmndFile)
        df₅ = DataFrame(csv₅)
        origin = df₅[!, 1]::Array{Int64,1}
        destination = df₅[!, 2]::Array{Int64,1}
        flows = df₅[!, 3:ncol(df₅)]::DataFrame
        for i in 1:nrow(df₅)
            rₒ = origin[i]
            for j in 1:(ncol(df₅)-2)
                if j == 1 r = origin[i]
                else r = length(N) + 1 end
                s = destination[i]
                m = j
                if r ∉ R Sᵣ[r] = [] end
                if r ∉ R append!(R, r) end
                append!(Sᵣ[r], s)
                Mᵣ[r] = m
                qᵣ[(r,s)] = flows[i,j]
                if j > 1
                    if r ∉ N
                        append!(N, r)
                        append!(A[rₒ], r), push!(A, [rₒ])
                        append!(Vᵢⱼ[rₒ], 1.0), push!(Vᵢⱼ, [1.0])
                        append!(dᵢⱼ[rₒ], 0.0), push!(dᵢⱼ, [0.0])
                        append!(tᵢⱼ[rₒ], 0.001), push!(tᵢⱼ, [0.001])
                        append!(αᵢⱼ[rₒ], 0.0), push!(αᵢⱼ, [0.0])
                        append!(βᵢⱼ[rₒ], 0.0), push!(βᵢⱼ, [0.0])
                        append!(ϕ[rₒ], 0), push!(ϕ, [0])
                    end
                end
            end
        end
    end

    # Returns cost of arc (i,j) for class m given arc flow x (k = A[i]⁻¹(j))
    function cᵢⱼ(i, k, m, x)
        #j = A[i][k]

        α = αᵢⱼ[i][k]
        β = βᵢⱼ[i][k]
        tₒ= tᵢⱼ[i][k]
        V = Vᵢⱼ[i][k]

        t = tₒ * (1 + α * (abs(x)/V) ^ β)
        d = dᵢⱼ[i][k]
        v = d/t
        if v == Inf v = 1.0e6 end

        c = 0.0
        for k in 1:length(η[m]) c += η[m][k] * v^(length(η[m]) - k) * t end
        return c
    end

    # Returns derivative of cost of arc (i,j) at arc flow x (k = A[i]⁻¹(j))
    function c′ᵢⱼ(i, k, m, x)
        #j = A[i][k]

        α = αᵢⱼ[i][k]
        β = βᵢⱼ[i][k]
        tₒ= tᵢⱼ[i][k]
        V = Vᵢⱼ[i][k]

        t = tₒ * (1 + α * (abs(x)/V) ^ β)
        d = dᵢⱼ[i][k]
        v = d/t
        if v == Inf v = 1.0e6 end

        t′ = tₒ * α * β * (abs(x) ^ (β - 1))/(V ^ β)
        if β == 0 t′ = 0.0 end
        if t′ == Inf t′ = 1.0e6 end

        c′ = 0.0
        for k in 1:length(η[m]) c′ += -(length(η[m]) - k - 1) * η[m][k] * v^(length(η[m]) - k) * t′ end
        return c′
    end

    # Returns cost for segment e given arc flows xₐ and arc costs c
    function cₑ(e, cₐ)
        c = 0.0
        for (n,i) in enumerate(e[1:end-1])
            j = e[n+1]
            k = findfirst(x -> (x == j), A[i])::Int64
            c += cₐ[i][k]
        end
        return c
    end

    # Returns flow on segment e given arc flows xₐ
    function fₑ(e, xₐ)
        f = Float64[]
        for (n,i) in enumerate(e[1:end-1])
            j = e[n+1]
            k = findfirst(x -> (x == j), A[i])::Int64
            append!(f, xₐ[i][k])
        end
        return minimum(f)
    end

    # Returns predecessor label L for every node i for least cost path from node r given arc costs cₐ
    function djk(cₐ, r)
        L = [if i==r r else -1 end for i in N]         # Predecessor label
        C = [if i==r 0.0 else Inf end for i in N]      # Cost label
        X = copy(N)                                    # Set of open nodes
        i = r
        deleteat!(X, i)
        while !isempty(X)
            for (k,j) in enumerate(A[i])
                c = C[i] + cₐ[i][k]
                if c < C[j] && j in X L[j], C[j] = i, c end
            end
            index = argmin([C[i] for i in X])
            i = X[index]
            deleteat!(X, index)
        end
        return L
    end

    # Returns tree rooted at r given predecessor label L
    function tree(L, r)
        T = Array{Int64,1}[[] for j in N]
        for j in N
            i = L[j]
            if i ≠ j && i ≠ -1 append!(T[i], j) end
        end
        return T
    end

    # Returns path between node r and s using predecessor label L
    function path(L, r, s)
        p = Int64[]
        i = s
        append!(p, i)
        while i ≠ r
            i = Int(L[i])
            append!(p, i)
        end
        reverse!(p)
        return p
    end

    # improved Traffic Assignment with Paired Alterantive Segments
    function iTAPAS(ϵ, θ, μ, 𝜈)
        report = Dict("TF" => Float64[], "TC" => Float64[], "RG" => Float64[], "WT" => Float64[])

        xʳₐ = Dict(r => [[0.0 for j in A[i]] for i in N] for r in R)                      # Stores origin-based arc flows
        xₐ  = [[sum([xʳₐ[r][i][k] for r in R]) for k in 1:length(A[i])] for i in N]       # Stores arc flows

        cₐ  = [[[cᵢⱼ(i, k, m, xₐ[i][k]) for k in 1:length(A[i])] for i in N] for m in M]     # Stores arc cost
        c′ₐ = [[[c′ᵢⱼ(i, k, m, xₐ[i][k]) for k in 1:length(A[i])] for i in N] for m in M]    # Stores derivative of arc cost
        πʳₐ = Dict(r => [[0.0 for j in A[i]] for i in N] for r in R)                         # Stores arc reduced cost

        rₚ = Int64[]                                                                        # Stores origin for PAS p
        P = Tuple{Array{Int64,1},Array{Int64,1}}[]                                          # Stores PAS

        Lᵣ = Dict(r => [if i==r r else -1 end for i in N] for r in R)                      # Stores origin-based least cost lables

        # Checks if arc a fails reduced cost optimal conditions for origin r
        function ispotential(a, r)
            i, j = a
            m = Mᵣ[r]
            k = findfirst(x -> (x == j), A[i])::Int64
            pᵣᵢ = path(Lᵣ[r], r, i)
            pᵣⱼ = path(Lᵣ[r], r, j)
            uʳᵢ = cₑ(pᵣᵢ, cₐ[m])
            uʳⱼ = cₑ(pᵣⱼ, cₐ[m])
            cʳᵢⱼ = cₐ[m][i][k]
            πʳᵢⱼ = uʳᵢ + cʳᵢⱼ - uʳⱼ
            if xʳₐ[r][i][k] > ϵ && πʳᵢⱼ > θ return (true, πʳᵢⱼ)
            else return (false, 0.0) end
        end

        # Checks if PAS p assosciated with origin rₒ can be eliminated
        function isbad(p, rₒ)
            e₁, e₂ = p
            m = Mᵣ[rₒ]
            c₁, c₂ = cₑ(e₁, cₐ[m]), cₑ(e₂, cₐ[m])
            f₁, f₂ = fₑ(e₁, xʳₐ[rₒ]), fₑ(e₂, xʳₐ[rₒ])
            if (f₁ < ϵ || f₂ < ϵ) && (c₁ ≠ c₂) return true
            else return false end
        end

        #= Shifts flows from higher cost segment to lower cost segment of PAS p
        on its assosciated origin rₒ, given cost difference is greater than λ=#
        function shift(p, rₒ, λ)
            e₁, e₂ = p
            m = Mᵣ[rₒ]

            c₁, c₂ = cₑ(e₁, cₐ[m]), cₑ(e₂, cₐ[m])
            if abs(c₂ - c₁) < λ return end

            c′₁, c′₂ = cₑ(e₁, c′ₐ[m]), cₑ(e₂, c′ₐ[m])
            f₁, f₂ = fₑ(e₁, xʳₐ[rₒ]), fₑ(e₂, xʳₐ[rₒ])
            Δ = (c₂ - c₁)/(c′₁ + c′₂)
            if isnan(Δ) δ = 0.0
            elseif Δ ≥ 0 δ = min(Δ, f₂)
            else δ = max(Δ, -f₁) end

            for (n,i) in enumerate(e₁[1:end-1])
                j = e₁[n+1]
                k = findfirst(x -> (x == j), A[i])::Int64
                xʳₐ[rₒ][i][k] += δ
                xₐ[i][k] += δ
                for m in M
                    cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k])
                    c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k])
                end
            end
            for (n,i) in enumerate(e₂[1:end-1])
                j = e₂[n+1]
                k = findfirst(x -> (x == j), A[i])::Int64
                xʳₐ[rₒ][i][k] -= δ
                xₐ[i][k] -= δ
                for m in M
                    cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k])
                    c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k])
                end
            end
        end

        # PAS identification for arc a wrt origin r using Maximum Cost Search
        function MCS(a, r)
            depth, maxdepth = 1, 2
            flag = true

            i, j = a
            e₁, e₂ = Int64[], Int64[]
            pᵣⱼ = path(Lᵣ[r], r, j)

            while flag
                # Intialize
                lₖ = [if k ∉ pᵣⱼ 0 elseif k ∈ a 1 else -1 end for k in N]
                L = [if k == j i else -1 end for k in N]

                # Iterate
                t = i
                while true
                    h = t

                    f = 0.0
                    for i in N
                        k = findfirst(x -> (x == h), A[i])
                        if k != nothing
                            x = xʳₐ[r][i][k]
                            m = Mᵣ[r]
                            c = cᵢⱼ(i, k, m, x)
                            if x > ϵ && c > f
                                f, t = c, i
                            end
                        end
                    end
                    L[h] = t
                    if lₖ[t] == -1      # PAS found
                        e₁ = path(Lᵣ[r], t, j)
                        e₂ = path(L, t, j)
                        shift((e₁, e₂), r, 0)
                        bool,_ = ispotential(a, r)
                        if !bool || depth == maxdepth flag = false
                        else depth += 1 end
                        break
                    elseif lₖ[t] == 1   # Cycle found
                        if depth == maxdepth flag = false
                        else
                            if h == t pₕₜ = Int64[]
                            else
                                pₕₜ = path(L, h, t)
                                append!(pₕₜ, h)
                                δ = fₑ(pₕₜ, xʳₐ[r])
                            end
                            for (n,i) in enumerate(pₕₜ[1:end-1])
                                k = findfirst(x -> (x == pₕₜ[n+1]), A[i])::Int64
                                xʳₐ[r][i][k] -= δ
                                xₐ[i][k] -= δ
                                for m in M
                                    cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k])
                                    c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k])
                                end
                            end
                            depth += 1
                        end
                        break
                    else                # Continue
                        lₖ[t] = 1
                    end
                end
            end
            p = (e₁, e₂)
            return p
        end

        ## Step 0: Intialization - AON assignment
        T =  Dates.format(now(), "HH:MM:SS:sss")
        tₒ = parse.(Int64, [T[1:2], T[4:5], T[7:8], T[10:12]])
        wt = 0.0
        iter = 0
        for r in R
            m = Mᵣ[r]
            Lᵣ[r] = djk(cₐ[m], r)
            for s in Sᵣ[r]
                qᵣₛ = qᵣ[r,s]
                pᵣₛ = path(Lᵣ[r], r, s)
                for (n,i) in enumerate(pᵣₛ[1:end-1])
                    j = pᵣₛ[n+1]
                    k = findfirst(x -> (x == j), A[i])
                    xʳₐ[r][i][k] += qᵣₛ
                    xₐ[i][k] += qᵣₛ
                    for m in M
                        cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k])
                        c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k])
                    end
                end
            end
        end

        if log == :on
            print("\n iter: iteration,  RG:Relative Gap,  TF:Total Flow,  TC: Total Cost,  WT: Wall Time (s)")
            print("\n iter  | logRG      | TF          | TC          | WT (s) ")
            print("\n ------|------------|-------------|-------------|--------")
        end

        # Iterate
        while true
            # Run Time calculation
            T =  Dates.format(now(), "HH:MM:SS:sss")
            tₙ = parse.(Int64, [T[1:2], T[4:5], T[7:8], T[10:12]])
            wt = sum((tₙ - tₒ) .* [3600, 60, 1, 1/1000])

            # Relative Gap calculation
            for r in R Lᵣ[r] = djk(cₐ[Mᵣ[r]], r) end
            num , den = 0.0, 0.0
            for r in R for s in Sᵣ[r] num += qᵣ[r,s] * cₑ(path(Lᵣ[r], r, s), cₐ[Mᵣ[r]]) end end
            for r in R for i in N for k in 1:length(A[i]) den += xʳₐ[r][i][k] * cₐ[Mᵣ[r]][i][k] end end end
            rg = 1 - num/den

            # Miscelaneous
            append!(report["RG"], log10(abs(rg)))
            append!(report["TF"], sum(sum.(xₐ)))
            append!(report["TC"], den)
            append!(report["WT"], wt)
            if log == :on
                if iter < 10 @printf("\n #%.0f    | %.3E | %.5E | %.5E | %.3f  ", iter, log10(abs(rg)), sum(sum.(xₐ)), den, wt)
                else @printf("\n #%.0f   | %.3E | %.5E | %.5E |%.3f ", iter, log10(abs(rg)), sum(sum.(xₐ)), den, wt) end
            end

            # Convergence Test
            if log10(abs(rg)) ≤ log10(tol) || iter ≥ maxIters || wt ≥ maxRunTime break end
            iter += 1

            ## Step 1
            for r in R
                m = Mᵣ[r]
                Lᵣ[r] = djk(cₐ[m], r)
                Tᵣ = tree(Lᵣ[r], r)
                ## Step 1.1: Indentify potential arcs
                for i in N
                    for (k,j) in enumerate(A[i])
                        a = (i, j)
                        bool, πʳᵢⱼ = false, 0.0
                        if j ∉ Tᵣ[i] bool, πʳᵢⱼ = ispotential(a, r) end
                        πʳₐ[r][i][k] = πʳᵢⱼ
                        ## Step 1.2: Flow shift for potential arc
                        if bool
                            p = MCS(a,r)
                            e₁, e₂ = p
                            rₒ = r
                            if !isempty(e₁) && !isempty(e₂) && p ∉ P push!(P, p), append!(rₚ, rₒ) end
                        end
                    end
                end
                # Step 1.3: Local shift
                for (k,p) in enumerate(P) shift(p, rₚ[k], rg/1000) end
            end
            ## Step 2
            for _ in 1:20
                for (k,p) in enumerate(P)
                    if isbad(p, rₚ[k]) deleteat!(P, k), deleteat!(rₚ, k)
                    else shift(p, rₚ[k], rg/1000) end
                end
            end
        end

        output = Dict("Flows" => xₐ, "Costs" => cₐ)
        println("\n Total run time: $wt; Total network flow: $(sum([xʳₐ[r][i][k] * ϕ[i][k] for r in R for i in N for k in 1:length(A[i])])); Total network cost: $(sum([xʳₐ[r][i][k] * cₐ[Mᵣ[r]][i][k]  * ϕ[i][k] for r in R for i in N for k in 1:length(A[i])]))")
        return output, report
    end
    build()
    return iTAPAS(1e-12, 1e-16, 0.5, 0.25)

end

traffic_assignment(networkName="Anaheim", tol=1e-12, maxIters=20, maxRunTime=600, log=:on)

# TODO: Combine the results in a post-processing step to render a multi-class solution - SemiDONE
# TODO: Upload networks in a repository and update repository name in the docstring
