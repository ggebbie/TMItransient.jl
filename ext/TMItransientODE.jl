module TMItransientODE

using OrdinaryDiffEq

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


"""
    function globalmean_rampresponse

calculate the global mean response to a ramp (linearly increasing value) in some region
"""
function globalmean_rampresponse(TMIversion, region, γ, L, B, τ)

    # assume evenly spaced (uniform) time spacing
    # Δτ = diff(τ)[1]
    c₀ = vec(zeros(γ)) # preallocate initial condition Field
    #q = vec(ones(γ)) # preallocate initial condition Field
    b = TMI.surfaceregion(TMIversion,region)
    q = B*vec(b)
    #qfunc = t -> q
    pfixed = [L,q,B.rowval]
    f(du,u,p,t) = constant_forcing!(du, u, pfixed) #avoid allocation
    #f(du,u,p,t) = mul!(du,  u, p[1] ) #avoid allocation
    func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array

    # make sure it starts at t=0 even if not saved there
    tspan = (0*first(τ),last(τ))
    #prob = ODEProblem(constant_forcing!, c₀, tspan, q) # Field type
    prob = ODEProblem(func, c₀, tspan) # Field type

    # possible algs:
    # QNDF, TRBDF2, FBDF, CVODE_BDF, lsoda, ImplicitEuler
    integrator = init(prob, QNDF())

    # better to grab input type somehow, instead of assuming Float64
    Dmean = Float64[] # [0.0]; # for time 0

    solfld = zeros(γ)
    for (u, t) in TimeChoiceIterator(integrator, τ)
        solfld.tracer[wet(solfld)] = u
        push!(Dmean,mean(solfld))
    end

    return Dmean
end


"""
    function globalmean_stepresponse_with_restoring

calculate the global mean response to an ocean restored to a step function
"""
function globalmean_stepresponse_with_restoring(TMIversion, region, γ, Lrestore, B, τ, τrestore)

    # assume evenly spaced (uniform) time spacing
    # Δτ = diff(τ)[1]
    c₀ = vec(zeros(γ)) # preallocate initial condition Field
    b = TMI.surfaceregion(TMIversion,region)
    θtarget = B*vec(b)

    # reset L in mixed layer or surface
    # for i in B.rowval
    #     L[i,i] = -1.0 / τrestore
    # end

    dutarget = (1.0/τrestore) * θtarget # set overriding and restoring boundary condition at right location.
    pfixed = (Lrestore, dutarget)

    jac_sparsity = ADTypes.jacobian_sparsity(
        (du,u) -> restored_forcing_test!(du, u , pfixed, 0.0), du, u, detector)

    #f(du,u,p,t) = restored_forcing!(du, u, Lrestore, dutarget) #avoid allocation
    #f(u,p,t) = restored_forcing(u, p, t) # take on some allocation
    #f(du,u,p,t) = mul!(du,  u, p[1] ) #avoid allocation
    f(du,u,p,t) = restored_forcing!(du, u, p, t) #avoid allocation

    #func = ODEFunction(f, jac = jacobian) #jac_prototype for sparse array
    func = ODEFunction(f, jac_prototype = float.(jac_sparsity)) #jac_prototype for sparse array
    println("func defined")
    
    # make sure it starts at t=0 even if not saved there
    tspan = (0*first(τ),last(τ))
    #prob = ODEProblem(constant_forcing!, c₀, tspan, q) # Field type
    prob = ODEProblem(func, c₀, tspan, pfixed) # Field type
    # prob = ODEProblem(func, c₀, tspan) # Field type

    # possible algs:
    # QNDF, TRBDF2, FBDF, CVODE_BDF, lsoda, ImplicitEuler
    integrator = init(prob, QNDF())
    #integrator = init(prob, TRBDF2())

    # better to grab input type somehow, instead of assuming Float64
    Dmean = Float64[] # [0.0]; # for time 0

    solfld = zeros(γ)
    for (u, t) in TimeChoiceIterator(integrator, τ)
        #solfld.tracer[wet(solfld)] = u
        #push!(Dmean,mean(solfld))
    end

    return Dmean
end

"""
    function globalmean_stepresponse_qndf

calculate the global mean response to "turning on" some region
"""
function globalmean_stepresponse_qndf(TMIversion,region,γ,L,B,τ)

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
    #integrator = init(prob,QNDF())
    integrator = init(prob,QNDF1())
    #integrator = init(prob,SSPSDIRK2())

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


end
