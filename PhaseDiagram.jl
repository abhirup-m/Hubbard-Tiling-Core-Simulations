using LinearAlgebra, Plots, ProgressMeter
include("Constants.jl")
include("RgFlow.jl")

function PoleFraction(
        size_BZ::Int64,
        J_val::Float64,
        W_val::Float64;
    )
    kondoJArray = momentumSpaceRG(size_BZ, J_val, W_val)
    _, dispersion = getDensityOfStates(tightBindDisp, size_BZ)
    fermiPoints = unique(getIsoEngCont(dispersion, 0.0))
    averageKondoScale = J_val / 2
    @assert averageKondoScale > RG_RELEVANCE_TOL
    kondoJArray .= ifelse.(abs.(kondoJArray) ./ averageKondoScale .> RG_RELEVANCE_TOL, kondoJArray, 0)
    polesFraction = count(p -> any(>(0), abs.(kondoJArray[p, fermiPoints])), fermiPoints)
    return polesFraction
end


function CriticalBathInt(
        size_BZ::Int64,
        omega_by_t::Float64,
        kondoJ::Float64,
        transitionWindow::Vector{Float64},
        bathIntSpacing::Float64;
        maxIter=100,
        loadData::Bool=false,
    )
    fracToIndex(f) = ifelse(f == 1, 1, ifelse(f > 0, 2, 3))

    savePathCrit = joinpath(SAVEDIR, "crit-$(size_BZ)-$(kondoJ).json")
    criticalBathIntData = Dict{String,Vector{Float64}}()
    if isfile(savePathCrit) && loadData
        merge!(criticalBathIntData, JSON3.read(read(savePathCrit, String), typeof(criticalBathIntData)))
        if string(bathIntSpacing) ∈ keys(criticalBathIntData)
            return criticalBathIntData[string(bathIntSpacing)]
        end
    end
    @assert bathIntSpacing > 0
    criticalBathInt = Float64[]
    @assert issorted(transitionWindow, rev=true)

    savePath = joinpath(SAVEDIR, "pf-$(size_BZ)-$(kondoJ).json")
    availableData = Dict{String,Float64}()
    if isfile(savePath)
        merge!(availableData, JSON3.read(read(savePath, String), typeof(availableData)))
    end
    for phaseBoundType in [(1, 2), (2, 3)]
        currentTransitionWindow = copy(transitionWindow)
        currentPoleFractions = [PoleFraction(size_BZ, omega_by_t, kondoJ, W_val; availableData=availableData, loadData=loadData) for W_val in currentTransitionWindow]
        currentPhaseIndices = map(fracToIndex, currentPoleFractions)
        numIter = 1
        while abs(currentTransitionWindow[1] - currentTransitionWindow[2]) > bathIntSpacing && numIter < maxIter
            updatedEdge = 0.5 * sum(currentTransitionWindow)
            newPoleFraction = PoleFraction(size_BZ, omega_by_t, kondoJ, updatedEdge; availableData=availableData, loadData=loadData)
            newPhaseIndex = fracToIndex(newPoleFraction)
            if newPhaseIndex == currentPhaseIndices[1] || newPhaseIndex == phaseBoundType[1]
                currentPhaseIndices[1] = newPhaseIndex
                currentTransitionWindow[1] = updatedEdge
            else
                currentPhaseIndices[2] = newPhaseIndex
                currentTransitionWindow[2] = updatedEdge
            end
            numIter += 1
        end
        push!(criticalBathInt, 0.5 * sum(currentTransitionWindow))
    end
    criticalBathIntData[string(bathIntSpacing)] = criticalBathInt

    open(savePath, "w") do file JSON3.write(file, availableData) end
    open(savePathCrit, "w") do file JSON3.write(file, criticalBathIntData) end

    return criticalBathInt
end

function PhaseDiagram(
        size_BZ::Int64,
        kondoJVals, 
        bathIntVals,
    )
    polesFractionMatrix = zeros(length(bathIntVals), length(kondoJVals))
    @showprogress for (x, J_val) in enumerate(kondoJVals)
        for (y, W_val) in enumerate(bathIntVals)
            polesFractionMatrix[y, x] = PoleFraction(size_BZ, J_val, W_val)
        end
    end
    p = heatmap(polesFractionMatrix, title=size_BZ)
    display(p)
end

PhaseDiagram(33, 0.05:0.01:0.3, 0:-0.01:-0.4)
