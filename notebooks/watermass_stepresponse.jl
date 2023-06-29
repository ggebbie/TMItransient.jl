### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2cdfb47e-15fc-11ee-096b-07dbb5eafd21
import Pkg; Pkg.activate(".")

# ╔═╡ 77f37139-fbd0-4a90-9273-0308a9b586fa
using TMItransient, TMI, Statistics, InteractiveUtils, Interpolations, PlutoUI, Plots, LinearAlgebra

# ╔═╡ 75aceba0-1532-4f44-8fd9-32fc89ac0482
md"""
# Water Mass Step Response Test 
This test demonstrates how to compute both the step response and impulse response to a boundary condition. The latter is computed by taking the derivative of the step response. 
"""

# ╔═╡ 8d59134d-a010-4c28-a5a4-c57100a352d6
begin
	TMIversion = "modern_90x45x33_GH10_GH12"
	#TMIversion = "modern_180x90x33_GH10_GH12"
	#TMIversion = "modern_90x45x33_unpub12"
	
	A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion,compute_lu=false);
end

# ╔═╡ 8732ad6e-a702-403b-b76d-a669eac8eba0
begin
	list = ("GLOBAL","ANT","SUBANT",
			"NATL","NPAC","TROP","ARC",
			"MED","ROSS","WED","LAB","GIN",
			"ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
			"TROPATL","TROPPAC","TROPIND")
	τ = 0:3
	@bind region Slider([l for l in list], show_value = true)
end

# ╔═╡ 299b8225-91c0-491f-8556-dd77439f9863
md"""
Calculate the step response to a surface condition in a particular region over the lags $\tau$. Similarly to `vintagetest.jl`, we generate a boundary condition and solve 

$$\frac{\partial c}{\partial t} = Lc$$

The function `globalmean_stepresponse` demonstrates how to solve the differential equation demonstrated here, but not save any output. It only computes statistics at the evenly-spaced times in $\tau$ for efficiency (it still is fairly slow). 

This function will output a vector, showing the mean percentage of interior waters sourced from the surface at each time in $\tau$. 
"""

# ╔═╡ c018635a-d366-4023-a1ef-a0a060c6b10b
#takes ~100-200 seconds to compute (harder to compute where the gradients are higher, so the first time step is very computationally expensive to compute.)
@time D̄ = globalmean_stepresponse(TMIversion,region,γ,L,B,τ)

# ╔═╡ d5f23f9e-39b0-49bf-a949-4a15d917c123
md"""
Next, we can calculate the impulse response by computing a simple derivative of the prior output. 
"""

# ╔═╡ 5b097a90-ee4f-4fd2-86bb-4fef0360965c
Ḡ,tḠ = globalmean_impulseresponse(TMIversion,region,γ,L,B,τ,alg=:centered)

# ╔═╡ Cell order:
# ╠═75aceba0-1532-4f44-8fd9-32fc89ac0482
# ╠═2cdfb47e-15fc-11ee-096b-07dbb5eafd21
# ╠═77f37139-fbd0-4a90-9273-0308a9b586fa
# ╠═8d59134d-a010-4c28-a5a4-c57100a352d6
# ╠═8732ad6e-a702-403b-b76d-a669eac8eba0
# ╠═299b8225-91c0-491f-8556-dd77439f9863
# ╠═c018635a-d366-4023-a1ef-a0a060c6b10b
# ╠═d5f23f9e-39b0-49bf-a949-4a15d917c123
# ╠═5b097a90-ee4f-4fd2-86bb-4fef0360965c
