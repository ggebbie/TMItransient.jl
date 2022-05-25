module TMItransient

using OrdinaryDiffEq, PreallocationTools, LinearAlgebra, NCDatasets, Interpolations

export readopt, ces_ncwrite, varying!, setupODE

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


end
