using TMItransient, TMI 
using Test

@testset "TMItransient.jl" begin
    # Write your tests here.
    
    TMIversion = "modern_90x45x33_GH10_GH12"
    #TMIversion = "modern_180x90x33_GH10_GH12"
    #TMIversion = "modern_90x45x33_unpub12"
    
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

    @testset "monotonicinterpolation" begin
        using Interpolations
        Δ,τ = read_stepresponse()
        itp = interpolate(τ, Δ, FritschCarlsonMonotonicInterpolation())
        itp = interpolate(τ, Δ, SteffenMonotonicInterpolation())
        #itp = interpolate(τ, Δ, FritschButlandInterpolation())

        #g = vintagedistribution(1850,2022,Δ,τ)

    end
    
    @testset "transientsimulation" begin

        using Interpolations, NaNMath, DifferentialEquations, LinearAlgebra, PreallocationTools

        latbox = [50,60]
        lonbox = [-50,0]
        d = surfacepatch(lonbox, latbox, γ) 
        dsfc = d.tracer[d.wet]

        #following make_initial_conditions.m
        c0 = B * dsfc 

        #Fixed euler timestep approximation
        c = c0
        Δt = 1e-3 #this becomes unstable if you go any lower
        T  = 1e-2
        Nt = T/Δt
        for tt = 1:Nt
            # forward Euler timestep
            c += L*c*Δt
            println("Σ c = ",sum(c))
        end
        gain_euler = sum(c .- c0)
        
        #Solving differential equation for fixed case 
        u0 = c0
        du = similar(u0)
        f(du,u,p,t) = mul!(du, L, u) 
        tspan = (0.0,T)
        func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
        prob = ODEProblem(func, u0, tspan)
        println("Solving fixed ODE")
        @time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,calck=false)
        println("ODE solved")

        #put sol into time x lon x lat x depth 
        sol_array = zeros((length(sol.t), 90,45,33))
        [sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]

        #stability check
        stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
        @test stable
        println("fixed bc stable: ", stable)

        #gain check - tracer concentration should increase 
        gain_ode = NaNMath.sum(sol_array[end, :, :, :].-sol_array[begin, :, :, :])
        println("Gain = ", gain_ode)
        @test gain_ode ≥ 0.0

        #compare forward euler timestep approx and solved ODE results 
        gain_error = abs(gain_ode - gain_euler)/(abs(gain_ode) + abs(gain_euler))
        @test gain_error < 0.1
        
        println("Gain percent error ",200gain_error,"%")

        #varying case stability check
        tsfc = [0, T]
        Csfc = zeros((2, length(dsfc)))
        Csfc[1, :] .= 1
        τ = 1/12
        li = LinearInterpolation(tsfc, 1:length(tsfc))
        LC = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl
        BF = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl 
        Cb = similar(Csfc[1,:])
        surface_ind = findall(x->x[3] ==1, γ.I)

        p = (Csfc,surface_ind,τ,L,B,li,LC,BF,Cb) #parameters
        f(du, u, p, t) = TMItransient.varying!(du, u, p, t)
        func = ODEFunction(f, jac_prototype=L)
        prob = ODEProblem(func, u0, tspan,p)
        println("Solving varying ODE")
        @time sol = solve(prob, QNDF(),abstol=1e-2,reltol=1e-2,saveat=tsfc)
        println("Varying ODE solved")
        
        #put sol into time x lon x lat x depth 
        sol_array = zeros((length(sol.t), 90,45,33))
        [sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]
        
        #stability check
        stable = true ? NaNMath.maximum(sol_array) < 1.000001 && NaNMath.minimum(sol_array) > -0.000001 : false
        @test stable
        println("Varying case stable: ", stable)
       
    end
end
