using Revise, TMItransient

Δ,τ = read_stepresponse()

g = vintage(1850,2022,Δ,τ)

# get some TTDs so that we can take difference of TTDs
