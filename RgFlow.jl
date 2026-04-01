using LinearAlgebra, Plots
include("Constants.jl")

######### HELPER FUNCTIONS #############
#
# helper functions for switching back and forth between the 
# 1D flattened representation (1 → N^2)
# and the 2D representation ((1 → N)×(1 → N))
function map1DTo2D(
        point::Int64,
        size_BZ::Int64
    )
    # Convert overall point to row and column values.
    # These serve as indices of kx and ky
    kx_index = (point - 1) % size_BZ + 1
    ky_index = div((point - 1), size_BZ) + 1

    # construct a 1D array of possible k values, and convert
    # kx_index and ky_index into values
    k_values = range(K_MIN, stop=K_MAX, length=size_BZ)
    return [k_values[kx_index], k_values[ky_index]]
end


function map1DTo2D(
        point::Vector{Int64},
        size_BZ::Int64
    )
    # same as above, but for multiple points. In this case,
    # two tuples are returned, for kx values and ky values.
    kx_index = (point .- 1) .% size_BZ .+ 1
    ky_index = div.((point .- 1), size_BZ) .+ 1
    k_values = range(K_MIN, stop=K_MAX, length=size_BZ)
    return [k_values[kx_index], k_values[ky_index]]
end


function map2DTo1D(
        kx_val::Float64,
        ky_val::Float64,
        size_BZ::Int64
    )
    # obtain the indices, given the values of kx and ky
    k_values = range(K_MIN, stop=K_MAX, length=size_BZ)
    kx_index, ky_index = argmin(abs.(k_values .- kx_val)), argmin(abs.(k_values .- ky_val))

    # covert the row and column indices into an overall flattened index
    return kx_index + (ky_index - 1) * size_BZ
end


function map2DTo1D(
        kx_val::Vector{Float64},
        ky_val::Vector{Float64},
        size_BZ::Int64
    )
    # same but for multiple (kx,ky) pairs.
    k_values = range(K_MIN, stop=K_MAX, length=size_BZ)
    points = Int64[]
    for (kx, ky) in zip(kx_val, ky_val)
        kx_index, ky_index = argmin(abs.(k_values .- kx)), argmin(abs.(k_values .- ky))
        push!(points, kx_index + (ky_index - 1) * size_BZ)
    end
    return points
end


# create flattened tight-binding dispersion
function tightBindDisp(kx_vals::Vector{Float64}, ky_vals::Vector{Float64})
    dispersion = -2 * HOP_T .* (cos.(kx_vals) + cos.(ky_vals))
    dispersion[abs.(dispersion) .< TOLERANCE] .= 0
    return dispersion
end
function tightBindDisp(kx_val::Float64, ky_val::Float64)
    dispersion = -2 * HOP_T .* (cos.(kx_val) + cos.(ky_val))
    dispersion[abs.(dispersion) .< TOLERANCE] .= 0
    return dispersion
end


function getDensityOfStates(
        dispersionFunc, 
        size_BZ::Int64,
    )
    kx_vals = repeat(range(K_MIN, stop=K_MAX, length=size_BZ), outer=size_BZ)
    ky_vals = repeat(range(K_MIN, stop=K_MAX, length=size_BZ), inner=size_BZ)

    dispersion = dispersionFunc(kx_vals, ky_vals)
    dispersion_xplus1 = dispersionFunc(circshift(kx_vals, 1), ky_vals)
    dispersion_xminus1 = dispersionFunc(circshift(kx_vals, -1), ky_vals)
    dispersion_yplus1 = dispersionFunc(ky_vals, circshift(ky_vals, size_BZ))
    dispersion_yminus1 = dispersionFunc(ky_vals, circshift(ky_vals, -size_BZ))
    dOfStates =
        4 ./
        sqrt.(
            (dispersion_xplus1 - dispersion_xminus1) .^ 2 +
            (dispersion_yplus1 - dispersion_yminus1) .^ 2
        )

    # van-Hoves might result in Nan, replace them with the largest finite value
    replace!(dOfStates, Inf => maximum(dOfStates[dOfStates.≠Inf]))

    # normalise the DOS according to the convention \int dE \rho(E) = N (where N is the total number of k-states)
    dOfStates /=
        sum([
            dos * abs(dispersion[i%size_BZ+1] - dispersion[i]) for
            (i, dos) in enumerate(dOfStates)
        ]) / size_BZ^2
    return dOfStates, dispersion
end


function getIsoEngCont(dispersion::Vector{Float64}, probeEnergy::Float64)
    # obtain the k-space points that have the specified energy. We count one point per row.

    contourPoints = Int64[]
    for (point, energy) in enumerate(dispersion)
        if abs(energy - probeEnergy) < TOLERANCE
            push!(contourPoints, point)
        end
    end
    return contourPoints
end


function getIsoEngCont(dispersion::Vector{Float64}, probeEnergies::Vector{Float64})
    # obtain the k-space points that have the specified energy. We count one point per row.

    contourPoints = Int64[]
    for (point, energy) in enumerate(dispersion)
        for probeEnergy in probeEnergies
            if abs(energy - probeEnergy) < TOLERANCE
                push!(contourPoints, point)
            end
        end
    end
    return contourPoints
end

function particleHoleTransf(
        point::Int64,
        size_BZ::Int64
    )
    # obtain the particle-hole transformed point k --> (k + π) % π.
    kx_val, ky_val = map1DTo2D(point, size_BZ)
    kx_new = kx_val <= 0.5 * (K_MAX + K_MIN) ? kx_val + 0.5 * (K_MAX - K_MIN) : kx_val - 0.5 * (K_MAX - K_MIN)
    ky_new = ky_val <= 0.5 * (K_MAX + K_MIN) ? ky_val + 0.5 * (K_MAX - K_MIN) : ky_val - 0.5 * (K_MAX - K_MIN)
    return map2DTo1D(kx_new, ky_new, size_BZ)
end

function particleHoleTransf(
        points::Vector{Int64},
        size_BZ::Int64
    )
    # obtain the particle-hole transformed point k --> (k + π) % π.
    kx_vals, ky_vals = map1DTo2D(points, size_BZ)
    kx_new = [kx <= 0.5 * (K_MAX + K_MIN) ? kx + 0.5 * (K_MAX - K_MIN) : kx - 0.5 * (K_MAX - K_MIN) for kx in kx_vals]
    ky_new = [ky <= 0.5 * (K_MAX + K_MIN) ? ky + 0.5 * (K_MAX - K_MIN) : ky - 0.5 * (K_MAX - K_MIN) for ky in ky_vals]
    return map2DTo1D(kx_new, ky_new, size_BZ)
end
######### END OF HELPER FUNCTIONS #############

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
            GMatrix[q, q] = densityOfStates[q] ./ (OMEGA_BY_t * HOP_T - energyCutoff / 2 + kondoJArray[q, q] / 4 + bathW / 2) 
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

function main()
    size_BZ = 77
    maps = []
    k1x_vals, k1y_vals = map1DTo2D(collect(1:size_BZ^2), size_BZ)

    # bare Kondo array for comparison
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
end

###### UNCOMMENT THE FOLLOWING LINE TO RUN THIS FILE ########
# main()
