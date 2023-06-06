using Revise, TMItransient

Δ,τ = read_stepresponse()

g = vintagedistribution(1850,2022,Δ,τ)

# next line not currently working
g = vintagedistribution(1850,2022,Δ,τ,tmodern=2022,interp="spline")
