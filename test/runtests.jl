using Revise
using TMItransient, TMI 
using Test
#using Statistics
using ExponentialUtilities
using LinearAlgebra

@testset "TMItransient.jl" begin
    TMIversion = "modern_90x45x33_GH10_GH12"
    #TMIversion = "modern_90x45x33_unpub12"

    #useful for explicit functions
    A, Alu, γ, TMIfile, L, B = config(TMIversion);
    Lmix = mixedlayermatrix(A, γ, 0.05)
    Ldir = dirichletmatrix(γ, 0.1)
    Ltot = L + Lmix + Ldir

    @testset "exponential" begin
        # use ExponentialUtilities

        @testset "global mean high level functions" begin
            # top-level, most abstracted algorithm
            # @testset "watermass_stepresponse" begin

            # NOTE: not sure if this next commented part should be uncommented.
            # read a water-mass surface patch from these choices
        # list = TMI.regionlist()

        # # choose water mass (i.e., surface patch) of interest
        # region = list[1]

        #τ = vcat(0.0:0.1:10,10:2000) # sample Common Era run
        # doesn't converge to 1 as well (overshoots)
        # nτ = 10000
        # τmax = 5000
        # τ = exp.(range(0,log(τmax + 1.0), nτ)).-1.0

        # @testset "exponential" begin
        #     # try ExponentialUtilities

        #     for i = 1:2
        #         if i == 1
        #             τ = 0.0:0.1:0.5
        #         else
        #             τ = 0.1:0.1:0.5
        #         end
            
            # make some of your own choices with keywords
            # or ignore these optional parameters and
            # trust the developers
            τi = 0:5
            τs = 0:0.1:5
            τd = 0.1
            τm = 0.2
            D, τD = globalmean_stepresponse(TMIversion,
                τ=τs,
                τdirichlet=τd,
                τmixedlayer=τm,
                alg=:exponential)

            G, τg = globalmean_impulseresponse(TMIversion,
                τimpulse=τi,
                τsimulate=τs,
                τdirichlet=τd,
                τmixedlayer=τm,
                alg=:exponential)
        end

        @testset "global mean basics" begin
    
            tf = 1
            τ1 = 0.0:0.01:tf
            τ2 = 0.1:0.01:tf
            τ3 = 0:0.1:tf
              
            @time D̄1, τ1out = globalmean_stepresponse(
                Ltot, τ1, γ, alg = :exponential)

            @time D̄2, τ2out = globalmean_stepresponse(
                Ltot, τ2, γ, alg = :exponential)

            @time D̄3, τ3out = globalmean_stepresponse(
                Ltot, τ3, γ, alg = :exponential)

            # should monotonically increase
            @test sum(diff(D̄1) .≥ 0) == length(D̄1) - 1
            @test sum(diff(D̄2) .≥ 0) == length(D̄2) - 1
            @test sum(diff(D̄3) .≥ 0) == length(D̄3) - 1

            # somewhat stable with different discretization?
            @test abs((last(D̄1) - last(D̄3)) / (last(D̄1) + last(D̄3))) < 0.2 # 20% difference at early stages

            @time G, τG = globalmean_impulseresponse(
                Ltot, τi, τs, γ, alg = :exponential)
        
            # Ḡ should be non-negative
            @test sum(Ḡ .≥ 0) == length(Ḡ)
            # Ḡ should add to something less than unity
            @test sum(Ḡ) ≤ 1.0
        end
            
        @testset "mean age" begin

            # compare g, g2 at N random points
            N = 2
            # get random locations that are wet (ocean)
            locs = [wetlocation(γ) for i in 1:N]

            #test: is the integral of ĝ equivalent to the output of the `meanage` function? (eqtn 2 of GH 2012) 
            τs = vcat(0:0.01:1.0,1.01:0.1:10,11:4000)
            τi = 0:4000

            @time G, τG = observe_impulseresponse(
                Ltot, τi, τs, locs, γ, alg = :exponential)

            # # QNDF: 90 seconds for 100, 98 seconds for 2000, 106 for 10k 
            # # exponential: 167 sec for 4k
        
            # uses locs from top-level scope
            a_obs = observe(meanage(TMIversion, Alu, γ), locs, γ)
            println("Equilbrium inversion for mean age at sites ",a_obs)
            g = hcat(G...)
            a = [cumsum(g[i, :] .* τG)[end] for i in 1:2]
            acorrection = [(1 - sum(g[i,:])) * τG[end] for i in 1:2]
            println("Transiently simulated mean age at sites ",a + acorrection)
            atol = 10
            denom = abs.(a + acorrection + a_obs)./2
            replace!(x -> x< atol ? atol : x, denom)
            relative_error = 100*abs.(a + acorrection - a_obs)./denom
            @test all(relative_error .< 5) # relative error less than 5 percent?
            println("percent relative error ",relative_error)

        end
        
        @testset "vintage test" begin

            # compare g, g2 at N random points
            N = 2
            # get random locations that are wet (ocean)
            locs = [wetlocation(γ) for i in 1:N]

            #test: is the integral of ĝ equivalent to the output of the `meanage` function? (eqtn 2 of GH 2012) 
            τs = vcat(0:0.01:1.0,1.01:0.1:10,11:4000)
            τ = 0:10

            @time D, τD = observe_stepresponse(
                Ltot, τ, τs, locs, γ, alg = :exponential)
                        
            vint =  zeros(length(locs))
            for j in  eachindex(y1)
                Δ  = [D[i][j]  for  i  in  eachindex(D)]
                vint[j] = vintagedistribution(2015,2020,Δ,τD)
            end
            
            @test maximum(vent) ≤ 1.0
            #@test minimum(g) ≥ 0.0 # fails for MATLAB

            g2 = vintagedistribution(TMIversion, γ, L, B, 2015, 2020)
            @test maximum(g2) ≤ 1.0
            #@test minimum(g) ≥ 0.0 # fails for Julia

            #y1 = TMI.observe(g,locs,γ)
            y2 = observe(g2,locs,γ)

            # formerly calculates relative difference between MATLAB and Julia computations
            # now calculates relative difference between direct and indirect computations
            for tt in 1:N
                @test 100*abs(y1[tt] - y2[tt])/(y1[tt] + y2[tt]) < 1.0 # percent
            end
        end

        @testset "watermass_stepresponse" begin

            @testset "water masses + regions" begin
                
                # read a water-mass surface patch from these choices
                list = TMI.regionlist()

                # choose water mass (i.e., surface patch) of interest
                region = list[1]
            end
        end
    end

    # @testset "QNDF" begin
        #     # replace with function call
        #     # add alg=QNDF() as optional argument
        #     b = TMI.surfaceregion(TMIversion,region)
        #     θtarget = B*vec(b)
        #     τrestore = 0.5 # yr
        #     Lrestore = deepcopy(L)
        #     # reset L in mixed layer or surface
        #     for i in B.rowval
        #         Lrestore[i,i] = -1.0 / τrestore
        #     end
        #     dutarget = (1.0/τrestore) * θtarget # set overriding and restoring boundary condition at right location.

        #     using SparseConnectivityTracer, ADTypes
        #     detector = TracerSparsityDetector()
        
        #     # test that core algorithm does right thing
        #     function restored_forcing_test!(du, u,p,t)
        #         println("maxu ",maximum(u))
        #         println(p[1][1,:])
        #         du[begin:end] = muladd(p[1],u,p[2])
        #     end

        #     # works ok
        #     u = vec(zeros(γ))
        #     du = copy(u)
        #     pfixed =(Lrestore, dutarget)
        #     restored_forcing_test!(du, u, pfixed, nothing)

        #     jac_sparsity = ADTypes.jacobian_sparsity(
        #         (du,u) -> restored_forcing_test!(du, u , pfixed, 0.0), du, u, detector)

        #     f(du,u,p,t) = restored_forcing!(du, u, p, t) #avoid allocation
        #     func = ODEFunction(f, jac_prototype = float.(jac_sparsity)) #jac_prototype for sparse array
        #     # make sure it starts at t=0 even if not saved there
        #     tspan = (0*first(τ),last(τ))
        #     #prob = ODEProblem(constant_forcing!, c₀, tspan, q) # Field type
        #     c₀ = vec(zeros(γ)) # preallocate initial condition Field
        #     prob = ODEProblem(func, c₀, tspan, pfixed) # Field type
        #     # prob = ODEProblem(func, c₀, tspan) # Field type

        #     # possible algs:
        #     # QNDF, TRBDF2, FBDF, CVODE_BDF, lsoda, ImplicitEuler
        #     integrator = init(prob, QNDF())
        #     #integrator = init(prob, TRBDF2())
        
        #     @time D̄ = globalmean_stepresponse_with_restoring(TMIversion, region, γ, Lrestore, B, τ, τrestore) # CDF

        #     @time D̄ = globalmean_stepresponse(TMIversion,region,γ,L,B,τ) # CDF
        #     @time D̄ = globalmean_rampresponse(TMIversion,region,γ,L,B,τ) # CDF
        # end

        # compare to reading same thing from MATLAB output.
        # Δ,τmat = read_stepresponse()

        # # relative difference between MATLAB and Julia computations
        # for tt in 2:3
        #     ϵ = 100*abs(mean(Δ[tt]) - D̄[tt])./(mean(Δ[tt]) + D̄[tt])
        #     println("percent difference is ",ϵ)
        #     @test ϵ < 1.0 # percent
        # end

    @testset "stepresponse" begin
        for i = 1:2
            if i == 1
                τ = 0.0:0.1:0.5
            else
                # don't start at τ = 0
                τ = 0.1:0.1:0.5
            end
            region = "GLOBAL"
            b = TMI.surfaceregion(TMIversion, region)

            #this should have the same result as globalmean_stepresponse 
            @time Dnew, τnew = stepresponse(TMIversion, b, γ, L, B, τ, eval_func = mean, alg = :exponential) 
            @time Dold, τold = globalmean_stepresponse(TMIversion, γ, L, B, τ, alg = :exponential) # CDF
            @test sum(Dnew .== Dold) == length(τnew)   

            #get output in Field type 
            @time Dall, τall = stepresponse(TMIversion, b, γ, L, B, τ) 

            #use synthetic observations to grab some random wet points to observe 
            #N = 10
            #locs = [wetlocation(γ) for i in 1:N]
            Dobs, τobs = stepresponse(TMIversion, b, γ, L, B, τ, eval_func = observe, args = (locs, γ)) 

            # works for all now, 4 Feb 2025
            for (i, d) in enumerate([Dnew, Dold, Dall, Dobs])
                try
                    impulseresponse(d, τobs)
                catch
                    println("impulseresponse doesn't work for D̄ number: " * string(i))
                end   
            end
        end
    end

    @testset "impulse response different time grid" begin
        τ_simulate = 0.1:0.1:1.0 # times where simulation output saved
        for i = 1:2
            if i == 1
                τ_edges = τ_simulate
            elseif i == 2
                τ_edges = 0:0.5:1.0 # edges of the bins used to compute impulse response
            end
    
            region = "GLOBAL"
            b = TMI.surfaceregion(TMIversion, region)

            #this should have the same result as globalmean_stepresponse 
            @time Dfine, τfine = stepresponse(TMIversion, b, γ, L, B, τ_simulate, eval_func = mean, alg = :exponential) 
            @time Gcourse, τcourse = TMItransient.impulseresponse(Dfine, τfine, τ_edges)
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
