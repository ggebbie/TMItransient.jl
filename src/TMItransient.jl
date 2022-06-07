module TMItransient

using OrdinaryDiffEq, PreallocationTools, LinearAlgebra, NCDatasets,
    MAT, TMI, Interpolations

export readopt, ces_ncwrite, varying!,
    read_stepresponse, datadir, srcdir, vintagedistribution

# Define these paths by hand so that we don't
# have to use DrWatson
pkgdir() = dirname(dirname(pathof(TMItransient)))
pkgdir(args...) = joinpath(pkgdir(), args...)

datadir() = joinpath(pkgdir(),"data")
datadir(args...) = joinpath(datadir(), args...)

srcdir() = joinpath(pkgdir(),"src")
srcdir(args...) = joinpath(srcdir(), args...)

"""
     read surface layer
"""
function readopt(filename,γ)
    nc = NCDataset(filename)
#    lat = nc["latitude"][:]
#    lon = nc["longitude"][:]
    time = nc["year"][:]
#    depth = nc["depth"][:]
    theta = nc["theta"][:, :, :, :]

    #flip to time descending order
    reverse!(time)
    time = convert(Vector{Int}, time)
    reverse!(theta, dims = 1)
    theta_permuted = zeros((size(theta)[1], size(γ.wet)[1], size(γ.wet)[2], size(γ.wet)[3]))
    permutedims!(theta_permuted, theta, [1,4,3,2])
    return time, theta_permuted 
end

"""
    function varying!(du, u, p, t)
    ODE function for varying boundary cond
    Sets up dc/dt = L*C + B*f to be solved 
# Arguments
- `du`: dc/dt (must have this name for DifferentialEquations.jl to work
- `u`: C, what we are solving for 
- `p`: parameters for diffeq - must hold specified vars  
- `t`: time we are solving for (automatically determined by DE.jl)
# Output
- `du`: numerical value of LC+Bf, vector of size 74064 for 4°
"""
function varying!(du, u, p, t)
    
    #load parameters 
    Csfc,surface_ind,τ,L,B,li,LC,BF,Cb = p
    println("time = ", t)
    
    #generate Cb - interpolated surface boundary condition  
    li_t = convert(Float64, li[t])
    Cb .= (ceil(li_t)-li_t).*Csfc[Int(floor(li_t)), :] .+ (li_t-floor(li_t)).*Csfc[Int(ceil(li_t)), :]
    
    #use PreallocationTools.jl to handle Dual type in u 
    LC = get_tmp(LC, first(u)*t) 
    BF = get_tmp(BF, first(u)*t) 

    #Figure out what u is at surface 
    u_sfc = @view u[surface_ind]

    #Inplace math to make faster 
    mul!(LC, L, u) 
    mul!(BF, B, -(u_sfc.-Cb)./τ) 
    @. du = LC + BF 
    nothing
end

"""
    function ces_ncwrite(γ,time,sol_array)
    Write .nc file output for commonerasim.jl 
# Arguments
- `γ`: 
- `time`: vector of time values 
- `sol_array`: solution array in form time x lat x lon x depth - must match γ + time 
# Output
- saves .nc file titled "ces_output.nc" in data array 
"""
function ces_ncwrite(γ,time,sol_array, filepath)
    file = filepath * "/ces_output.nc"
    ds = NCDataset(file,"c")

    #define dimensions 
    defDim(ds,"lon", size(γ.lon)[1])
    defDim(ds,"lat",size(γ.lat)[1])
    defDim(ds,"depth",size(γ.depth)[1])
    defDim(ds,"time",size(time)[1])

    #write theta output variable 
    v = defVar(ds,"theta",Float64, ("time","lon","lat","depth"))
    v[:,:,:,:] = sol_array

    #write dimensions as variables (this might not be kosher...) 
    vlon = defVar(ds,"lon",Float64, ("lon",))
    vlon[:] = γ.lon
    vlat = defVar(ds,"lat",Float64,("lat",))
    vlat[:] = γ.lat
    vtime = defVar(ds,"time",Float64,("time",))
    vtime[:] = time
    vdepth = defVar(ds,"depth",Float64,("depth",))
    vdepth[:] = γ.depth

    v.attrib["title"] = "output of commonerasim.jl" 
    v.attrib["units"] = "potential temperature anomaly"
    close(ds)
end

""" 
    function gooddata
    a useful one-liner
"""
goodtime = x -> (typeof(x) <: Number && !isnan(x))

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

    Δ = Vector{Field{Float64}}(undef,ngood)
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
        Δ[tt] = Field(tracer,γ,:Δ,"step response","dimensionless")
        
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
    function vintagedistribution(t₀,tf,Δ,τ)

    percentage of water in the modern ocean
    from a vintage defined by the calendar year interval
    t₀ [cal yr CE] => tf [cal yr CE]

# Arguments
- `t₀`: Starting calendar year of vintage, e.g., 1450 CE
- `tf`: final calendar year of vintage, e.g., 1850 CE
- `Δ::Vector{Field}`: Step function response
- `τ = Vector{Float64}`: time lags associated with step response
- `tmodern=2022`: modern calendar year

# Output
- `g::Field`: 3D distribution of vintage contribution

# Warning
- should be a way to make Δ argument more general (more types)
"""
function vintagedistribution(t₀,tf,Δ,τ,tmodern=2022)

    τ₀ = tmodern - t₀ # transfer starting cal year to equivalent lag
    τf = tmodern - tf # end year

    interp_linear = LinearInterpolation(τ, Δ)
    #g = interp_linear(τ₀) - interp_linear(τf)

    if isinf(t₀)
        # we know that CDF(τ=Inf) = 1.
        return g = ones(Δ[1].γ) - interp_linear(τf)
    else
        return g = interp_linear(τ₀) - interp_linear(τf)
    end

end


    
end
