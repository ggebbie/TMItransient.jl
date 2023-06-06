using TMItransient, TMI 
using Test
using Statistics

@testset "TMItransient.jl" begin
    # Write your tests here.
    
    TMIversion = "modern_90x45x33_GH10_GH12"
    #TMIversion = "modern_180x90x33_GH10_GH12"
    #TMIversion = "modern_90x45x33_unpub12"
    
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion,compute_lu=false);

    @testset "vintage test" begin

        using Interpolations
        Δ,τ = read_stepresponse()
        g = vintagedistribution(2015,2020,Δ,τ)

        @test maximum(g) ≤ 1.0
        #@test minimum(g) ≥ 0.0 # fails for MATLAB


        g2 = vintagedistribution(TMIversion,γ,L,B,2015,2020)
        @test maximum(g2) ≤ 1.0
        #@test minimum(g) ≥ 0.0 # fails for Julia

        # compare g, g2 at N random points

        N = 2
        # get random locations that are wet (ocean)
        locs = Vector{Tuple{Float64,Float64,Float64}}(undef,N)
        [locs[i] = wetlocation(γ) for i in eachindex(locs)]

        # get weighted interpolation indices
        wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
        [wis[i] = interpindex(locs[i],γ) for i in 1:N]

        y1 = observe(g,wis,γ)
        y2 = observe(g2,wis,γ)

        # relative difference between MATLAB and Julia computations
        for tt in 1:N
            @test 100*abs(y1[tt] - y2[tt])/(y1[tt] + y2[tt]) < 1.0 # percent
        end

    end

    # @testset "monotonicinterpolation" begin
    #     using Interpolations
    #     Δ,τ = read_stepresponse()

    #     d2 = Vector{Vector{Float64}}(undef,2)
    #     d2 = Vector{Float64}(undef,2)
    #     t2 = Vector{Float64}(undef,2)
    #     d2[1] = [1,2]
    #     d2[2] = [3,4]
    #     d2[1] = 1
    #     d2[2] = 3
    #     t2[1] = 1
    #     t2[2] = 2
        
    #     itp = interpolate(t2, d2, FritschCarlsonMonotonicInterpolation())
    #     itp = interpolate(τ, Δ)
    #     itp = interpolate(τ, Δ, FritschCarlsonMonotonicInterpolation())
    #     itp = interpolate(τ, Δ, SteffenMonotonicInterpolation())
    #     #itp = interpolate(τ, Δ, FritschButlandInterpolation())

    #     #g = vintagedistribution(1850,2022,Δ,τ)
    # end

    @testset "watermass_stepresponse" begin
        #using Revise, TMI
        using LinearAlgebra
        #using OrdinaryDiffEq
        #using Interpolations
        #using PyPlot
        #using NaNMath
        #using PythonPlot

        #TMIversion = "modern_90x45x33_GH10_GH12"
        #A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

        # read a water-mass surface patch from these choices
        list = ("GLOBAL","ANT","SUBANT",
                "NATL","NPAC","TROP","ARC",
                "MED","ROSS","WED","LAB","GIN",
                "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
                "TROPATL","TROPPAC","TROPIND")

        # choose water mass (i.e., surface patch) of interest
        region = list[1]
        #tspan = (0.0, 5.0)
        τ = 0:6
        # replace with function call
        # add alg=QNDF() as optional argument

        @time D̄ = globalmean_stepresponse(TMIversion,region,γ,L,B,τ) # CDF

        # should monotonically increase
        @test sum(diff(D̄) .≥ 0) == length(D̄) - 1

        Ḡ,tḠ = globalmean_impulseresponse(TMIversion,region,γ,L,B,τ,alg=:centered)
        #Ḡ2,tḠ2 = globalmean_impulseresponse(TMIversion,region,γ,L,B,τ,alg=:leapfrog)
        #Ḡ3,tḠ3 = globalmean_impulseresponse(D̄,τ,alg=:leapfrog)
        #jldsave("Gbar_TMI_90x45x33_GH10_GH12.jld2";Ḡ,tḠ)
        
        # Ḡ should be non-negative
        @test sum(Ḡ .≥ 0) == length(Ḡ)

        # Ḡ should add to something less than unity
        @test sum(Ḡ) ≤ 1.0

        # compare to reading same thing from MATLAB output.
        Δ,τmat = read_stepresponse()

        # relative difference between MATLAB and Julia computations
        for tt in 2:3
            @test 100*abs(mean(Δ[tt]) - D̄[tt])./(mean(Δ[tt]) + D̄[tt]) < 1.0 # percent
        end
        
    end

    # @testset "transientsimulation" begin

    #     using Interpolations, NaNMath, DifferentialEquations, LinearAlgebra, PreallocationTools, Sundials

    #     latbox = [50,60]
    #     lonbox = [-50,0]
    #     d = surfacepatch(lonbox, latbox, γ) 
    #     dsfc = d.tracer[d.wet]

    #     #following make_initial_conditions.m
    #     c0 = B * dsfc 

    #     #Fixed euler timestep approximation
    #     c = c0
    #     Δt = 1e-3 #this becomes unstable if you go any lower
    #     T  = 1e-2
    #     Nt = T/Δt
    #     for tt = 1:Nt
    #         # forward Euler timestep
    #         c += L*c*Δt
    #         println("Σ c = ",sum(c))
    #     end
    #     gain_euler = sum(c .- c0)
        
    #     #Solving differential equation for fixed case 
    #     u0 = c0
    #     du = similar(u0)
    #     f(du,u,p,t) = mul!(du, L, u) 
    #     tspan = (0.0,T)
    #     func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array
    #     @testset "fixed ODE" begin

    #         prob = ODEProblem(func, u0, tspan)
    #         println("Solving fixed ODE")
    #         @time sol = solve(prob,CVODE_BDF(linear_solver=:GMRES),abstol = 1e-4,reltol=1e-4,calck=false)
    #         println("ODE solved")
        
    #         #put sol into time x lon x lat x depth 
    #         sol_array = zeros((length(sol.t), 90,45,33))
    #         [sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]

    #         stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
    #         println("fixed bc stable: ", stable)
    
    #         #gain check - tracer concentration should increase 
    #         gain_ode = NaNMath.sum(sol_array[end, :, :, :].-sol_array[begin, :, :, :])
    #         println("Gain = ", gain_ode)
    #         @test gain_ode ≥ 0.0

    #         #compare forward euler timestep approx and solved ODE results 
    #         gain_error = abs(gain_ode - gain_euler)/(abs(gain_ode) + abs(gain_euler))
    #         @test gain_error < 0.1
        
    #         println("Gain percent error ",200gain_error,"%")

    #         #varying case stability check
    #         tsfc = [0, T]
    #         Csfc = zeros((2, length(dsfc)))
    #         Csfc[1, :] .= 1
    #         τ = 1/12
    #         li = LinearInterpolation(tsfc, 1:length(tsfc))
    #     #LC = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl
    #     #BF = DiffEqBase.dualcache(similar(u0)) #for PreallocationTools.jl 
    #         LC = dualcache(similar(u0)) #for PreallocationTools.jl
    #         BF = dualcache(similar(u0)) #for PreallocationTools.jl 
    #         Cb = similar(Csfc[1,:])
    #         surface_ind = findall(x->x[3] ==1, γ.I)

    #         p = (Csfc,surface_ind,τ,L,B,li,LC,BF,Cb) #parameters
    #         f(du, u, p, t) = TMItransient.varying!(du, u, p, t)
    #         func = ODEFunction(f, jac_prototype=L)
    #         prob = ODEProblem(func, u0, tspan,p)
    #         println("Solving varying ODE")
    #         @time sol = solve(prob, QNDF(),abstol=1e-2,reltol=1e-2,saveat=tsfc)
    #         println("Varying ODE solved")
            
    #         #put sol into time x lon x lat x depth 
    #         sol_array = zeros((length(sol.t), 90,45,33))
    #         [sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]
            
    #         #stability check
    #         stable = true ? NaNMath.maximum(sol_array) < 1.000001 && NaNMath.minimum(sol_array) > -0.000001 : false
    #         @test stable
    #         println("Varying case stable: ", stable)
    #     end       
    #end
end
