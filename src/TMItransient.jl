module TMItransient

using OrdinaryDiffEq, PreallocationTools, LinearAlgebra

export ces_ncwrite, varying!

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
function ces_ncwrite(γ,time,sol_array)
    file = pkgdatadir() * "/ces_output.nc"
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
