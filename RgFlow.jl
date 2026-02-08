using LinearAlgebra, Plots
include("Constants.jl")
include("Helpers.jl")


function highLowSeparation(
        dispersionArray,
        energyCutoff,
        proceedFlags,
        size_BZ
    )

    # get the k-points that will be decoupled at this step, by getting the isoenergetic contour at the cutoff energy.
    cutoffPoints = unique(getIsoEngCont(dispersionArray, energyCutoff))
    cutoffHolePoints = particleHoleTransf(cutoffPoints, size_BZ)

    # these cutoff points will no longer participate in the RG flow, so disable their flags
    proceedFlags[[cutoffPoints; cutoffHolePoints], :] .= 0
    proceedFlags[:, [cutoffPoints; cutoffHolePoints]] .= 0

    # get the k-space points that need to be tracked for renormalisation, by getting the states 
    # below the cutoff energy. We only take points within the lower left quadrant, because the
    # other quadrant is obtained through symmetry relations.
    innerIndices = [
                    point for (point, energy) in enumerate(dispersionArray) if
                    abs(energy) < (abs(energyCutoff) - TOLERANCE)
                    && any(proceedFlags[:, point])
                   ]
    return innerIndices, cutoffPoints, cutoffHolePoints, proceedFlags
end


function momentumSpaceRG(
        size_BZ::Int64,
        kondoJ::Float64,
        bathW::Float64,
    )
    kvals = map1DTo2D.(1:size_BZ^2, size_BZ)
    kxVals = first.(kvals)
    kyVals = last.(kvals)
    omega_by_t = OMEGA_BY_t

    # ensure that [0, \pi] has odd number of states, so 
    # that the nodal point is well-defined.
    @assert (size_BZ - 5) % 4 == 0 "Size of Brillouin zone must be of the form N = 4n+5, n=0,1,2..., so that all the nodes and antinodes are well-defined."

    densityOfStates, dispersionArray = getDensityOfStates(tightBindDisp, size_BZ)

    kx_pos_arr = [kx for kx in range(K_MIN, K_MAX, length=size_BZ) if kx >= 0]
    cutOffEnergies = sort(-tightBindDisp(kx_pos_arr, 0 .* kx_pos_arr), rev=true)

    # Kondo coupling must be stored in a 2D matrix. The two dimensions store the 
    # incoming and outgoing momentum indices.For example, J[i][j] reveals the value of J 
    # for the momentum pair (k_i, k_j).
    k1x_vals, k1y_vals = map1DTo2D(collect(1:size_BZ^2), size_BZ)
    kondoJArray = 0.5 * kondoJ .* (cos.(k1x_vals' .- k1x_vals) .+ cos.(k1y_vals' .- k1y_vals))

    initSigns = sign.(kondoJArray)

    # define flags to track whether the RG flow for a particular J_{k1, k2} needs to be stopped 
    # (perhaps because it has gone to zero, or its denominator has gone to zero). These flags are
    # initialised to one, which means that by default, the RG can proceed for all the momenta.
    proceedFlags = fill(true, size_BZ^2, size_BZ^2)

    initDeltaSign = repeat([1.], size_BZ^2, size_BZ^2)

    WMatrix = 0.5 .* (cos.(kxVals' .- kxVals) .+ cos.(kyVals' .- kyVals))

    # Run the RG flow starting from the maximum energy, down to the penultimate energy (ΔE), in steps of ΔE
    for (stepIndex, energyCutoff) in enumerate(cutOffEnergies[1:end-1])
        deltaEnergy = abs(cutOffEnergies[stepIndex+1] - cutOffEnergies[stepIndex])

        # if there are no enabled flags (i.e., all are zero), stop the RG flow
        if !any(proceedFlags)
            break
        end

        # innerIndices is the set of points that belong in IR
        # cutoffPoints, cutoffHolePoints are the set of points that belong in the shell being decoupled
        # proceedFlags decides whether certain IR states have reached their fixed points
        innerIndices, cutoffPoints, cutoffHolePoints, proceedFlags = highLowSeparation(dispersionArray, energyCutoff, proceedFlags, size_BZ)
        GMatrix = 0 .* kondoJArray

        # calculate G_q = ρ_q/denominator_q for every state q on the shell
        for q in cutoffPoints
            GMatrix[q, q] = densityOfStates[q] ./ (omega_by_t * HOP_T - energyCutoff / 2 + kondoJArray[q, q] / 4 + bathW / 2) 
        end

        # calculates sum_q(J_{q,qbar} * G_q)
        traceGprime = sum([GMatrix[q, q] * kondoJArray[q, qbar] for (q, qbar) in zip(cutoffPoints, cutoffHolePoints)])

        # calculates ΔJ = -ΔE * (J_{k,q} * G_q * J_{q,k'} - 4 * sum_q(J_{q,qbar} * G_q) * W_{k,k'}
        delta = -abs(deltaEnergy) * (kondoJArray[innerIndices, cutoffPoints] * GMatrix[cutoffPoints, cutoffPoints] * kondoJArray[cutoffPoints, innerIndices] .- 4 * bathW * traceGprime .* WMatrix[innerIndices, innerIndices])

        # if the first step, observe the sign of the renormalisation,
        # otherwie check if the sign has changed compared to first step,
        # if it has then set the renormalisation to zero because we
        # must have gone through a fixed pont.
        if step == 1
            initDeltaSign = sign.(delta)
        else
            initDeltaSign[innerIndices, innerIndices][sign.(delta) .* initDeltaSign[innerIndices, innerIndices] .< 0] .= 0.
            delta[initDeltaSign[innerIndices, innerIndices] .== 0] .= 0
        end

        # apply the renormalisation, and check whether any coupling 
        # J_{ki,kj} has changed sign. If it has, set its value to 0
        # and set its proceedFlag to false because we don't need to
        # track it anymore.
        kondoJArray[innerIndices, innerIndices] .+= delta
        for (i, j) in Iterators.product(innerIndices, innerIndices)
            if kondoJArray[i, j] * initSigns[i, j] ≤ 0
                kondoJArray[i, j] = 0
                proceedFlags[i, j] = false
            end
        end
    end
    return kondoJArray
end

size_BZ = 77
maps = []
k1x_vals, k1y_vals = map1DTo2D(collect(1:size_BZ^2), size_BZ)

# bare kOndo array for comparison
initArray = 0.5 * 0.1 .* (cos.(k1x_vals' .- k1x_vals) .+ cos.(k1y_vals' .- k1y_vals))
averageKondoScale = sum(abs.(initArray)) / length(initArray)

densityOfStates, dispersionArray = getDensityOfStates(tightBindDisp, size_BZ)
FS = getIsoEngCont(dispersionArray, 0.0)
for W in -1 .* [0., 14.04, 14.5, 14.99] ./ 77
    @time kondoJ = momentumSpaceRG(size_BZ, 0.1, W)
    kondoJ[abs.(kondoJ) ./ averageKondoScale .< 1e-2] .= 0
    kondoDiagMatrix = zeros(size_BZ^2)
    kondoDiagMatrix[FS] = [sum(kondoJ[k, FS].^2)^0.5 / sum(initArray[k, FS].^2)^0.5 for k in FS]
    push!(maps, heatmap(reshape(kondoDiagMatrix, size_BZ, size_BZ), title="W/J=$(W/0.1)"))
end
display(plot(maps...))
