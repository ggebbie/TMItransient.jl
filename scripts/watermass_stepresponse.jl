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
#using LinearAlgebra
#using OrdinaryDiffEq
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

# do numerical analysis
b = TMI.surfaceregion(TMIversion,region,γ)

# preallocate initial condition Field
c₀ = zeros(γ)

# update d with the boundary condition b
#setboundarycondition!(c,b)

c₀.tracer[TMI.wet(c₀)] = B* b.tracer[b.wet]
#initial conditions are the surface patch = 1, propagated down to the bottom of the mixed layer, which we get from B
#c0 = B * dsfc 

u0 = c₀
du = zeros(γ)
tspan = (0.0, 50.0)

#Solving differential equation
#NOTE: for DifferentialEquations.jl to work, follow naming conventions in docs

f(du,u,p,t) = mul!(du, L, u) #this avoids allocating a new array for each iteration

#Solve diff eq
operator = DiffEqArrayOperator(L)
#isconstant(operator) # not currently working
func = ODEFunction(f, jac_prototype = L) #jac_prototype for sparse array 
prob = ODEProblem(func, u0, tspan)
#prob = ODEProblem(func, u0, tspan,p)
println("Solving ode")
#solve using QNDF alg - tested against other alg and works fastest 
@time sol = solve(prob,QNDF(),abstol = 1e-4,reltol=1e-4,saveat =tspan[1]:tspan[2])
println("ode solved")

#put sol into time x lon x lat x depth 
sol_array = zeros((length(sol.t),size(γ.wet)[1],size(γ.wet)[2],size(γ.wet)[3]))
[sol_array[i,:,:,:] = vec2fld(sol.u[i],γ.I) for i ∈ 1:length(sol.t)]

#stability check
stable = true ? NaNMath.maximum(sol_array) < 1.000001  && NaNMath.minimum(sol_array) > -0.000001 : false
println("stable: " *string(stable))

#____PLOTTING____
#longitudinal plots
#lon_index = 85
#dyeplot(γ.lat, γ.depth, sol_array[25, lon_index, :, :]', 0:0.05:1.05)

lon_section = 330; # only works if exact
lims = 0:0.05:3.0
snapshot = TMI.Field(sol_array[25,:,:,:],γ)
label = "Some quantity [units]"
sectionplot(snapshot,lon_section,lims,titlelabel=label)

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330 # only works if exact
lims = 0:5:100
tlabel = region * " water-mass fraction [%]"
sectionplot(100g,lon_section,lims,titlelabel = tlabel)
