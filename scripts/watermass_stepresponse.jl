#=
 Example 8: Propagate initial conditions through time according to L-matrix 
 Steps: (a) Define a surface patch (d) of concentration 1 to propagate
        (b) Use a fixed euler timestep approximation to estimate how L
            will propagate d
        (c) Solve dc/dt = Lu for the fixed boundary condition case 
 Varying boundary condition: Solve for tracer concentration given an initial
 global surface tracer concentration of 1 that linearly decreases for 50 years
 until global surface tracer concentration is 0. Functionality can be turned on
 by toggling bc from "fixed" to "varying" 
 Solves dc/dt = Lu + Bf 
=#

using Revise, TMI
using LinearAlgebra
#using OrdinaryDiffEq
#using Statistics
using Plots
using TMItransient
#using Interpolations
#using PyPlot
#using NaNMath
#using PythonPlot

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

# read a water-mass surface patch from these choices
list = ("GLOBAL","ANT","SUBANT",
            "NATL","NPAC","TROP","ARC",
            "MED","ROSS","WED","LAB","GIN",
            "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
            "TROPATL","TROPPAC","TROPIND")

# choose water mass (i.e., surface patch) of interest
region = list[1]
tspan = (0.0, 4000.0)

# replace with function call
# add alg=QNDF() as optional argument
Gmean = globalmean_stepresponse(TMIversion,region,γ,L,B,tspan)

# lando: neither plotlyjs nor gr working

## end of tested work

# do numerical analysis
b = TMI.surfaceregion(TMIversion,region,γ)

# preallocate initial condition Field
c₀ = zeros(γ)

# update d with the boundary condition b
#setboundarycondition!(c,b)

#c₀ = B* b.tracer[b.wet]
c₀ = B* vec(b)
#initial conditions are the surface patch = 1, propagated down to the bottom of the mixed layer, which we get from B
#c0 = B * dsfc 

#u0 = c₀
#du = zeros(γ)
tspan = (0.0, 1.0)

#Solving differential equation
f(du,u,p,t) = mul!(du, L, u) #avoid allocation

#Solve diff eq
#operator = DiffEqArrayOperator(L); # too big of an initial shock = unstable
#isconstant(operator) # not currently working
#prob = ODEProblem(operator, u0.tracer[TMI.wet(u0)], tspan) # too big of an initial shock

func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
prob = ODEProblem(func, c₀, tspan) # Field type
#prob = ODEProblem(func, u0.tracer[TMI.wet(u0)], tspan) # Vector type
#prob = ODEProblem(func, c₀.tracer[wet(c₀)], tspan) # 
# QNDF alg tested and fastest

# possible algs:
# QNDF, TRBDF2, FBDF, CVODE_BDF, lsoda, ImplicitEuler
@time sol = solve(prob,TRBDF2(),abstol = 1e-4,reltol=1e-4,saveat =tspan[2])
println("ode solved")

## only really care about global-response to global mean.

# put sol into Field
solfld = zeros(γ)
τ = 1:tspan[2]
Gmean = zeros(length(τ))
for (ii,tt) in enumerate(τ)
    solfld.tracer[wet(solfld)] = sol(tt)
    Gmean[ii] = mean(solfld)
end

## get global mean CDF
gr()
plot(τ,Gmean)

## take differences of global mean CDF to get G_mean



#timeseries = [sol.u[t][50000] for t = eachindex(sol.u)]

#put sol into time x lon x lat x depth 
# sol_array = zeros((length(sol.t),size(γ.wet)[1],size(γ.wet)[2],size(γ.wet)[3]))
# [sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]

# #stability check
# #stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
# #println("stable: " *string(stable))

# #____PLOTTING____
# #longitudinal plots
# #lon_index = 85
# #dyeplot(γ.lat, γ.depth, sol_array[25, lon_index, :, :]', 0:0.05:1.05)

# lon_section = 330; # only works if exact
# lims = 0:0.05:3.0
# snapshot = TMI.Field(sol_array[25,:,:,:],γ)
# label = "Some quantity [units]"
# sectionplot(snapshot,lon_section,lims,titlelabel=label)

# # plot a section at 330 east longitude (i.e., 30 west)
# lon_section = 330 # only works if exact
# lims = 0:5:100
# tlabel = region * " water-mass fraction [%]"
# sectionplot(100g,lon_section,lims,titlelabel = tlabel)
