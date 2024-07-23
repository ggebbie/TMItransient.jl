module TMItransient

using OrdinaryDiffEq
using PreallocationTools
using LinearAlgebra
using NCDatasets
using TMI
using Interpolations
using Statistics
#using MAT
#using NaNMath

export readopt, ces_ncwrite, varying!,
    setupODE, setupODE_nojac, s_array,
    vintagedistribution, agedistribution,
    EvolvingField,
    globalmean_stepresponse,
    globalmean_impulseresponse,
    stepresponse,  deltaresponse
#  read_stepresponse, 
#  deltaresponse, taudeltaresponse,
#  stability_check,

"""
    struct EvolvingField

    This structure permits the grid to be 
    automatically passed to functions with
    the evolving tracer field.

    This structure assumes the Tracer type to be 
    three-dimensional with an additional vector dimension for time.

    tracer::Vector{Array{T,3}}
    γ::Grid
    name::Symbol
    longname::String
    units::String
"""
struct EvolvingField{T}
    tracer::Vector{Array{T,3}}
    γ::Grid
    name::Symbol
    longname::String
    units::String
end

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
    #li = LinearInterpolation(tsfc, 1:length(tsfc))
    li = linear_interpolation(tsfc, 1:length(tsfc))

    #Instantiate arrays that the diffeq solver will reallocate
    LC = PreallocationTools.dualcache(similar(u0)) #for PreallocationTools.jl
    BF = PreallocationTools.dualcache(similar(u0)) #for PreallocationTools.jl 
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
    li = linear_interpolation(tsfc, 1:length(tsfc))

    #Instantiate arrays that the diffeq solver will reallocate
    LC = PreallocationTools.dualcache(similar(u0)) #for PreallocationTools.jl
    BF = PreallocationTools.dualcache(similar(u0)) #for PreallocationTools.jl 
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
        itp = linear_interpolation(τ, Δ)
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

# Try to compute without using MATLAB file
function vintagedistribution(TMIversion,γ::TMI.Grid,L,B,t₀,tf; tmodern= 2023)

    τ₀ = tmodern - t₀ # transfer starting cal year to equivalent lag
    τf = tmodern - tf # end year

    #  Δτ = diff(τ)[1]

    # vintages defined relative to GLOBAL surface
    b = TMI.surfaceregion(TMIversion,"GLOBAL")
    
    #c₀ = zeros(γ) # preallocate initial condition Field
    c₀ = B* vec(b)
    f(du,u,p,t) = mul!(du, L, u) #avoid allocation
    func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array
    tspan = (0.0,τ₀)
    prob = ODEProblem(func, c₀, tspan) # Field type

    # possible algs:
    # QNDF, TRBDF2, FBDF, CVODE_BDF, lsoda, ImplicitEuler
    u = solve(prob,QNDF(),saveat=(τf,τ₀))

    g = zeros(γ)
    g.tracer[wet(g)] = u[2] - u[1]

    return g

    # ### old version
    
    # # get interpolation object
    # if interp == "linear"
    #     itp = linear_interpolation(τ, Δ)
    # elseif interp == "spline"
    #     println("warning: not implemented yet for type Field")
    #     itp = interpolate(τ, Δ, FritschCarlsonMonotonicInterpolation())
    # end
    
    # if isinf(t₀)
    #     # we know that CDF(τ=Inf) = 1.
    #     # return g = ones(Δ[1].γ) - interp_linear(τf)
    #     return g = ones(Δ[1].γ) - itp(τf)
    # else
    #     #return g = interp_linear(τ₀) - interp_linear(τf)
    #     return g = itp(τ₀) - itp(τf)
    # end
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

"""
    function deltaresponse

    Take CDF and turn it into PDF
"""
function deltaresponse(Δ,τΔ)

    # hard coded annual resolution, easy because don't have to normalize
    τmax = maximum(τΔ)
    tgedge = 0:floor(τmax)
    tg = 0.5:floor(τmax)

    #    interp_linear = linear_interpolation(τΔ, Δ)
    # Δhires = interp_linear(tgedge)

    itp = interpolate(τΔ, Δ, FritschCarlsonMonotonicInterpolation())
    Δhires = itp.(tgedge)
    g = diff(Δhires)
        
    return tg, g
end

"""
    function stepresponse

    calculate the response to "turning on" some region
    can compute some statistics on output by providing a function to f 

    # Arguments
    - TMIversion
    - b: BoundaryCondition
    - γ
    - L
    - B
    - τ: evenly spaced vector
    - f: some function f(u) where u is a vector of all wet points 
"""
function stepresponse(TMIversion, b, γ, L, B, τ; eval_func = return_self, args = [])
    # assume evenly spaced (uniform) time spacing
    Δτ = diff(τ)[1]
    #b = TMI.surfaceregion(TMIversion,region,γ)
    c₀ = zeros(γ) # preallocate initial condition Field
    c₀ = B * vec(b)
    f(du,u,p,t) = mul!(du, L, u) #avoid allocation
    func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array
    tspan = (first(τ), last(τ))
    prob = ODEProblem(func, c₀, tspan) #Field type

    # possible algs:
    # QNDF, TRBDF2, FBDF, CVODE_BDF, lsoda, ImplicitEuler
    integrator = init(prob,QNDF())
    
    #assumes `f` returns one output!
    #how should I handle the fact that there can be no args
    output = isempty(args) ? Vector{first(Base.return_types(eval_func, (Field{Float64},)))}(undef, length(τ)) : Vector{first(Base.return_types(eval_func, (Field{Float64}, typeof.(args)...,)))}(undef, length(τ))
    
    solfld = zeros(γ) #initialize solution Field 
    
    for (idx, (u, t)) in enumerate(TimeChoiceIterator(integrator, τ))
        solfld.tracer[wet(solfld)] = u
        output[idx] = isempty(args) ? eval_func(solfld) : eval_func(solfld, args...)
    end
    return output
        
end

return_self(x) = x 

"""
function globalmean_stepresponse

    calculate the global mean response to "turning on" some region
"""
function globalmean_stepresponse(TMIversion,region,γ,L,B,τ)

    # assume evenly spaced (uniform) time spacing
    # Δτ = diff(τ)[1]
    b = TMI.surfaceregion(TMIversion,region)
    c₀ = zeros(γ) # preallocate initial condition Field
    c₀ = B* vec(b)
    f(du,u,p,t) = mul!(du, L, u) #avoid allocation
    func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array
    # make sure it starts at t=0 even if not saved there
    tspan = (0*first(τ),last(τ))
    prob = ODEProblem(func, c₀, tspan) # Field type

    # possible algs:
    # QNDF, TRBDF2, FBDF, CVODE_BDF, lsoda, ImplicitEuler
    integrator = init(prob,QNDF())

    # better to grab input type somehow, instead of assuming Float64
    Dmean = Float64[] # [0.0]; # for time 0

    solfld = zeros(γ)
    for (u,t) in TimeChoiceIterator(integrator,τ)
        solfld.tracer[wet(solfld)] = u
        push!(Dmean,mean(solfld))
    end

    # Philosophy: would prefer to not mess with output.
    #set first element to zero if lag is zero
    #if iszero(τ[1])
    #    Dmean[1] = 0.0
    #end
    
    return Dmean
end

#globalmean_impulseresponse(TMIversion,region,γ,L,B,τ) = (diff(globalmean_stepresponse(TMIversion,region,γ,L,B,τ)),(τ[1:end-1]+τ[2:end])./2)
"""
    function globalmean_impulseresponse

    Ḡ: satisfied ∫₀^∞ Ḡ(τ) dτ = 1

    `alg`: centered or leapfrog

    Does leapfrog satisfy normalization?
"""
globalmean_impulseresponse(TMIversion,region,γ,L,B,τ;alg=:centered) = globalmean_impulseresponse(globalmean_stepresponse(TMIversion,region,γ,L,B,τ),τ,alg=alg)

"""
    function globalmean_impulseresponse

    based on a globalmean_stepresponse D̄, compute the impulse response
    can be done via centered or leapfrog difference 
"""
function globalmean_impulseresponse(D̄,τ;alg=:centered)
    if alg == :centered
        ihi = 2:length(D̄)
        ilo = 1:(length(D̄)-1)
    elseif alg == :leapfrog
        ihi = 3:length(D̄)
        ilo = 1:(length(D̄)-2)
    end
    Δτ = τ[ihi]-τ[ilo]
    Ḡ = (D̄[ihi] - D̄[ilo])./Δτ
    τ = (τ[ihi] + τ[ilo])./2
    return Ḡ,τ
end

end
