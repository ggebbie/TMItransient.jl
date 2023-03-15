module TMItransient

using OrdinaryDiffEq
using PreallocationTools
using LinearAlgebra
using NCDatasets
using MAT
using TMI
using NaNMath
using Interpolations

export readopt, ces_ncwrite, varying!,
    setupODE,setupODE_nojac, s_array, stability_check,
    read_stepresponse, vintagedistribution,
    deltaresponse, taudeltaresponse, agedistribution

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
    #println("changes")
    #load parameters 
    Csfc,surface_ind,τ,L,B,li,LC,BF,Cb = p
    #println("time = ", t)
    
    #generate Cb - interpolated surface boundary condition  
    li_t = convert(Float16, li[t])
    csfc_floor = @view Csfc[Int(floor(li_t)), :]
    csfc_ceil = @view Csfc[Int(ceil(li_t)), :]
    frac = ceil(li_t) - li_t
    Cb .= frac.*csfc_floor .+ (1-frac).*csfc_ceil
    
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
    function setupODE(γ, u0,tsfc, dsfc,bc,L, t_int)
    Setup ODEFunction for LC+Bf 
# Arguments
- `γ`: from TMI 
- `u0`: initial condition
- `tsfc`: time that surface conditions are specfiied at 
- `dsfc`: 
- `bc`:4D array that we will take the surface behavior from to use as forcing
- `L`:from TMI 
- `t_int`: an integer in 1:length(tsfc) with how many timesteps we want to take
# Output
- `func`: ODEFunction object
"""

function setupODE(γ, u0,tsfc,dsfc,bc,L,B,t_int) 
    #Timespan that diffeq solver will solve for, must be within tsfc 
    tspan = (tsfc[begin], tsfc[t_int])

    #Get surface boundary conditions from Theta_anom_OPT 
    Csfc = zeros((length(tsfc), length(dsfc)))
    [Csfc[i,:] = bc[i,:,:,1][γ.wet[:,:,1]] for i ∈ 1:length(tsfc)]

    #more parameters for diffeq solver 
    τ = 1 / 12 #monthly restoring timescale
    li = LinearInterpolation(tsfc, 1:length(tsfc))

    #Instantiate arrays that the diffeq solver will reallocate
    LC = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl
    BF = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl 
    Cb = similar(Csfc[1,:])
    surface_ind = findall(x-> x[3] == 1, γ.I) #Find which points in γ.I are on the surface
    #setup ODEproblem and return 
    
    f(du, u, p, t) = varying!(du, u, p, t) #diffeq function to solve
    jacobian(u,p,t) = L 
    func = ODEFunction(f, jac = jacobian, jac_prototype = L) #jac_prototype for sparse array
    p = (Csfc, surface_ind, τ, L, B, li, LC, BF, Cb) 
    return func, p, tspan
end

"""
    function setupODE(γ, u0,tsfc, dsfc,bc,L, t_int)
    Setup ODEFunction for LC+Bf with no Jacobian 
"""
function setupODE_nojac(γ, u0,tsfc,dsfc,bc,L,B,t_int) 
    #Timespan that diffeq solver will solve for, must be within tsfc 
    tspan = (tsfc[begin], tsfc[t_int])

    #Get surface boundary conditions from Theta_anom_OPT 
    Csfc = zeros((length(tsfc), length(dsfc)))
    [Csfc[i,:] = bc[i,:,:,1][γ.wet[:,:,1]] for i ∈ 1:length(tsfc)]

    #more parameters for diffeq solver 
    τ = 1 / 12 #monthly restoring timescale
    li = LinearInterpolation(tsfc, 1:length(tsfc))

    #Instantiate arrays that the diffeq solver will reallocate
    LC = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl
    BF = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl 
    Cb = similar(Csfc[1,:])
    surface_ind = findall(x-> x[3] == 1, γ.I) #Find which points in γ.I are on the surface
    #setup ODEproblem and return 
    
    f(du, u, p, t) = varying!(du, u, p, t) #diffeq function to solve
    #jacobian(u,p,t) = L 
    func = ODEFunction(f)#, jac = jacobian, jac_prototype = L) #jac_prototype for sparse array
    p = (Csfc, surface_ind, τ, L, B, li, LC, BF, Cb) 
    return func, p, tspan
end

"""
    function  s_array(sol, γ)
    Converts from DE.jl output to time x lat x lon x depth 
"""
function s_array(sol, γ)
    sol_array = zeros((length(sol.t),size(γ.wet)[1],size(γ.wet)[2],size(γ.wet)[3]))
    [sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]
    return sol_array
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
    function vintagedistribution(t₀,tf,Δ,τ,tmodern=2022,interp="linear")

    percentage of water in the modern ocean
    from a vintage defined by the calendar year interval
    t₀ [cal yr CE] => tf [cal yr CE]

# Arguments
- `t₀`: Starting calendar year of vintage, e.g., 1450 CE
- `tf`: final calendar year of vintage, e.g., 1850 CE
- `Δ::Vector{Field}`: Step function response
- `τ = Vector{Float64}`: time lags associated with step response
- `tmodern=2022`: modern calendar year
- `interp=linear`: or can be "spline"
# Output
- `g::Field`: 3D distribution of vintage contribution

# Warning
- should be a way to make Δ argument more general (more types)
"""
function vintagedistribution(t₀,tf,Δ,τ;tmodern=2023,interp="linear")

    τ₀ = tmodern - t₀ # transfer starting cal year to equivalent lag
    τf = tmodern - tf # end year

    # get interpolation object
    if interp == "linear"
        itp = LinearInterpolation(τ, Δ)
    elseif interp == "spline"
        println("warning: not implemented yet for type Field")
        itp = interpolate(τ, Δ, FritschCarlsonMonotonicInterpolation())
    end
    
    if isinf(t₀)
        # we know that CDF(τ=Inf) = 1.
        # return g = ones(Δ[1].γ) - interp_linear(τf)
        return g = ones(Δ[1].γ) - itp(τf)
    else
        #return g = interp_linear(τ₀) - interp_linear(τf)
        return g = itp(τ₀) - itp(τf)
    end
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
    function agedistribution

    age distribution refers to distribution over lags, not space
    sometimes called a transit time distribution
"""
function agedistribution(loc)

    Δloc,τ = stepresponse(loc)
    tg, g = deltaresponse(Δloc,τ)
    
    return g
end

#     Δ,τ = read_stepresponse()
#     ncdf = length(τ)
#     npdf = ncdf - 1

#     # get weighted interpolation indices
#     wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,1)
#     #[wis[xx] = interpindex(loc[xx],Δ[xx].γ) for xx in 1]
#     wis[1] = interpindex(loc,Δ[1].γ)
    
#     Δloc = Vector{Float64}(undef,ncdf)
#     for tt in 1:ncdf
#         Δloc[tt] = observe(Δ[tt],wis,Δ[tt].γ)[1] # kludge to convert to scalar
#     end
#     println(typeof(Δloc))
#     tg, g = deltaresponse(Δloc,τ)
    
#     return g
    
# end

"""
    function deltaresponse

    Take CDF and turn it into PDF
"""
function deltaresponse(Δ,τΔ)

    # hard coded annual resolution, easy because don't have to normalize
    τmax = maximum(τΔ)
    tgedge = 0:floor(τmax)
    tg = 0.5:floor(τmax)

    #    interp_linear = LinearInterpolation(τΔ, Δ)
    # Δhires = interp_linear(tgedge)

    itp = interpolate(τΔ, Δ, FritschCarlsonMonotonicInterpolation())
    Δhires = itp.(tgedge)
    g = diff(Δhires)
        
    return tg, g
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
    stable = true ? NaNMath.maximum(sol_array) < NaNMath.maximum(Csfc) && NaNMath.minimum(sol_array) > NaNMath.minimum(Csfc) : false
    println("stable: " *string(stable))
end


end
