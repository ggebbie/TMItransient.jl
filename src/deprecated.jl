"""
    stability_check(sol_array, Csfc) 
    checks stability of ODE output 
# Arguments
- `sol_array`: solution array 
- `Csfc`: boundary condition 

# Output
- prints true or false 
"""
function stability_check(sol_array, Csfc) 
    #stable = true ? NaNMath.maximum(sol_array) < NaNMath.maximum(Csfc) && NaNMath.minimum(sol_array) > NaNMath.minimum(Csfc) : false
    println("stable: " *string(stable))
end

"""
    function read_stepresponse()

    read previously computed MATLAB output using shell script

    Warning: hard-coded file names from a 4° x 4° TMI run
"""
function read_stepresponse()

    TMIversion = "modern_90x45x33_GH10_GH12"

    ncfile = download_ncfile(TMIversion)
    γ = Grid(ncfile)

    # next load pre-computed step function response
    stepfile =  datadir("ttdsummary_global_13july2011.mat")
    !isfile(stepfile) && download_stepresponse()

    matfile = download_matfile(TMIversion)
    Izyx = cartesianindex(matfile)

    matobj = matopen(stepfile)
    y = read(matobj,"Y") # CDF of step function response
    
    τ = vec(read(matobj,"T")) # time lag

    # eliminate empty times
    ngood = count(goodtime,sum(y,dims = 2))

    Δ = Vector{TMI.Field}(undef,ngood)
    τ = τ[1:ngood]
    
    for tt in 1:ngood

        if iszero(τ[tt])
            # set to zero so that mixed-layer goes from zero to one
            # in the first year
            tracer = tracerinit(zeros(sum(γ.wet)), Izyx, γ.wet)

        else
            # initialize a 3D tracer array where zyx format
            # is transferred to xyz format
            tracer = tracerinit(y[tt,:], Izyx, γ.wet)
        end
        
        # construct a Field type
        Δ[tt] = TMI.Field(tracer,γ,:Δ,"step response","dimensionless")
        
    end
    
    return Δ, τ
end

function download_stepresponse()
    shellscript = srcdir("read_stepresponse.sh")
    run(`sh $shellscript`)

    # move the file to datadir
    println(joinpath(pwd(),"ttdsummary_global_13july2011.mat"))
    !isfile(joinpath(pwd(),"ttdsummary_global_13july2011.mat"))

    !isdir(datadir()) && mkdir(datadir())
    mv(joinpath(pwd(),"ttdsummary_global_13july2011.mat"),
       datadir("ttdsummary_global_13july2011.mat"))

end


"""
    function stepresponse(loc)

    response to a Heaviside function at a `loc`
"""
function stepresponse(loc)

    Δ,τ = read_stepresponse()
    ncdf = length(τ)
    npdf = ncdf - 1

    # get weighted interpolation indices
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,1)
    #[wis[xx] = interpindex(loc[xx],Δ[xx].γ) for xx in 1]
    wis[1] = interpindex(loc,Δ[1].γ)
    
    Δloc = Vector{Float64}(undef,ncdf)
    for tt in 1:ncdf
        Δloc[tt] = observe(Δ[tt],wis,Δ[tt].γ)[1] # kludge to convert to scalar
    end
    
    return Δloc,τ    
end


"""
    function taudeltaresponse

    Take CDF and turn it into PDF, get lag timescale

    Seems inefficient, reading file for time lag alone
"""
function taudeltaresponse()

    Δ,τ = read_stepresponse()

    # hard coded annual resolution, easy because don't have to normalize
    τmax = maximum(τ)
    tgedge = 0:floor(τmax)
    tg = 0.5:floor(τmax)
    return tg
end
