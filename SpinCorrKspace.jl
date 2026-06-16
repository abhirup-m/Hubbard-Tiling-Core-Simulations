using LinearAlgebra, Plots, ProgressMeter, Fermions
include("Constants.jl")
include("RgFlow.jl")


function KondoModel(
        dispersion::Vector{Float64},
        kondoJ::Matrix{Float64},
    )

    #### SITE INDEXING SCHEME
    # d↑    d↓    k1↑    k1↓    k2↑    k2↓
    # 1     2     3      4      5      6
    # ⟹    k_n↑ = 1 + 2 * n
    numBathSites = length(dispersion)
    hamiltonian = Tuple{String, Vector{Int64}, Float64}[]

    # kinetic energy
    for site in 1:numBathSites
        push!(hamiltonian, ("n",  [1 + 2 * site], dispersion[site])) # up spin
        push!(hamiltonian, ("n",  [2 + 2 * site], dispersion[site])) # down spin
    end

    # kondo terms
    for indices in Iterators.product(1:numBathSites, 1:numBathSites)
        up1, up2 = 2 .* indices .+ 1
        down1, down2 = (up1, up2) .+ 1
        push!(hamiltonian, ("n+-",  [1, up1, up2], kondoJ[indices...] / 4)) # n_{d up, n_{0 up}
        push!(hamiltonian, ("n+-",  [1, down1, down2], -kondoJ[indices...] / 4)) # n_{d up, n_{0 down}
        push!(hamiltonian, ("n+-",  [2, up1, up2], -kondoJ[indices...] / 4)) # n_{d down, n_{0 up}
        push!(hamiltonian, ("n+-",  [2, down1, down2], kondoJ[indices...] / 4)) # n_{d down, n_{0 down}
        push!(hamiltonian, ("+-+-",  [1, 2, down1, up2], kondoJ[indices...] / 2)) # S_d^+ S_0^-
        push!(hamiltonian, ("+-+-",  [2, 1, up1, down2], kondoJ[indices...] / 2)) # S_d^- S_0^+
    end

    return hamiltonian
end


function IterDiagKspace(
        hamiltDetails::Dict,
        maxSize::Int64,
        momentumPoints::Vector{Int64},
        pivot::Int64,
        correlationFuncDict::Dict{String, Any};
    )
    momentumPoints = sort(momentumPoints, by=p -> sum((map1DTo2D(p, hamiltDetails["size_BZ"]) .- map1DTo2D(pivot, hamiltDetails["size_BZ"])).^2))
    correlationDefDict = Dict{String, Vector{Tuple{String, Vector{Int64}, Float64}}}()
    for (name, func) in correlationFuncDict
        correlationDefDict[name] = func(3, 3)
    end
    hamiltonian = KondoModel(
                             hamiltDetails["dispersion"][momentumPoints],
                             hamiltDetails["kondoJArray"][momentumPoints, momentumPoints];
                            )
    indexPartitions = [10]
    while indexPartitions[end] < 2 + 2 * length(momentumPoints)
        push!(indexPartitions, indexPartitions[end] + 2)
    end
    hamiltonianFamily = MinceHamiltonian(hamiltonian, indexPartitions)

    results = IterDiag(
                      hamiltonianFamily, 
                      maxSize;
                      symmetries=Char['N', 'S'],
                      correlationDefDict=correlationDefDict,
                      maxMaxSize=2 * maxSize,
                     )

    return Dict(k => results[k] for k in keys(correlationDefDict))
end


function SpinCorrKspace(
        size_BZ::Int,
        bathIntVal::Float64,
        maxSize::Int,
    )
    J_val = 0.1
    correlationFuncDict = Dict{String, Any}("SPM-dk" => (i, j) -> [("+-+-", [1, 2, i+1, j], 1.0), ("+-+-", [2, 1, i, j+1], 1.0)])

    kondoJArray = momentumSpaceRG(size_BZ, J_val, bathIntVal)

    # identify points on Fermi surface. This is independent
    # of (J, W) and simply looks at the zero energy surface
    # of the non-interacting dispersion
    _, dispersion = getDensityOfStates(tightBindDisp, size_BZ)
    fermiPoints = unique(getIsoEngCont(dispersion, 0.0))

    # if certain J_{k,k'} matrix elements have fallen below
    # a threshold at the fixed point, set them to zero,
    # the threshold RG_RELEVANCE_TOL is set in Constants.jl
    averageKondoScale = J_val / 2
    @assert averageKondoScale > RG_RELEVANCE_TOL
    kondoJArray .= ifelse.(abs.(kondoJArray) ./ averageKondoScale .> RG_RELEVANCE_TOL, kondoJArray, 0)

    hamiltDetails = Dict()
    hamiltDetails["size_BZ"] = size_BZ
    hamiltDetails["kondoJArray"] = kondoJArray
    hamiltDetails["dispersion"] = dispersion
    corr = zeros(size_BZ^2)
    corr .= NaN
    pivotPoints = filter(p -> map1DTo2D(p, size_BZ)[1] ≥ π/2 && map1DTo2D(p, size_BZ)[2] ≥ 0, fermiPoints)
    @showprogress desc="overall: " Threads.@threads for pivot in pivotPoints
        results = IterDiagKspace(hamiltDetails, maxSize, fermiPoints, pivot, correlationFuncDict)
        corr[pivot] = results["SPM-dk"]
    end
    nodes = [(π/2, π/2), (π/2, -π/2), (-π/2, -π/2), (-π/2, π/2), (π, π), (π, -π), (-π, -π), (-π, π)]
    distances = [sum([sum((map1DTo2D(p, size_BZ) .- node).^2) for node in nodes])^2 for p in 1:size_BZ^2]
    for pivot in pivotPoints
        equidistant = sortperm(abs.(distances .- distances[pivot]))
        equidistant = filter(p -> dispersion[p] == dispersion[pivot], equidistant)[1:8]
        corr[equidistant] .= corr[pivot]
    end
    p = plot()
    heatmap!(p, reshape(corr, (size_BZ, size_BZ)))
    savefig(p, "SpinCorr.pdf")
    display(p)
end

### Run this if you DON'T have time
SpinCorrKspace(13, -0.15, 600)

### Run this if you have time
# SpinCorrKspace(33, -0.21, 600)
