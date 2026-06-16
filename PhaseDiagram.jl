using LinearAlgebra, Plots, ProgressMeter
include("Constants.jl")
include("RgFlow.jl")

# Given a set of parameters (J, W), simulates the RG
# equations for those parameters, and counts the number
# of points N_p on Fermi surface that satisfy |J_{k,k}| > 0
function CountPoles(
        size_BZ::Int64,
        J_val::Float64,
        W_val::Float64;
    )

    # obtain fixed point Kondo matrix for given (J, W)
    kondoJArray = momentumSpaceRG(size_BZ, J_val, W_val)

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

    # count the number of points {k} on Fermi surface that satisfy |J_kk| > 0.
    # this is taken to be the criteria for the said point remaining gapless.
    poles = count(p -> any(>(0), abs.(kondoJArray[p, fermiPoints])), fermiPoints)

    return poles
end


# Given a range of values for J and W, counts the 
# number of Kondo relevant Fermi surface points 
# for every value and plots it as a colourmap on J X W
function PhaseDiagram(
        size_BZ::Int64,
        kondoJVals, 
        bathIntVals,
    )
    # matrix in the space of J_i x W_j. Stores the number of
    # gapless points on Fermi surface (see `CountPoles` for
    # details) at any given parameter point (J_i, W_j)
    polesMatrix = zeros(length(bathIntVals), length(kondoJVals))

    # This can be parallelised for large gains. 
    # Uncomment the line below and comment out 
    # the line below it to enable multithreading.
    # Julia must be started with somthing like 
    # `julia -t 5` to leverage this, where 5 is 
    # the number of threads.
    
    # @showprogress Threads.@threads for ((x, J_val), (y, W_val)) in Iterators.product(enumerate(kondoJVals), enumerate(bathIntVals)) |> collect
    @showprogress for ((x, J_val), (y, W_val)) in Iterators.product(enumerate(kondoJVals), enumerate(bathIntVals))
        polesMatrix[y, x] = CountPoles(size_BZ, J_val, W_val)
    end
    p = heatmap(kondoJVals, abs.(bathIntVals), polesMatrix, title="\$N = $size_BZ\$", xlabel="\$J\$", ylabel="\$W\$")
    display(p)
    savefig(p, "phaseDiagram.pdf")
end

#### if you dont' have time
PhaseDiagram(13, 0.05:0.01:0.3, 0:-0.02:-0.4)

#### if you DO have time
# PhaseDiagram(33, 0.05:0.01:0.3, 0:-0.02:-0.4)
