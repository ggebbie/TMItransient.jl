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

# ╔═╡ 50b8ca96-1c49-4690-a0ae-5f8cacaea429
import Pkg; Pkg.activate(".")

# ╔═╡ a088c589-906d-4d18-b9d9-44dab6c7fbcd
using TMItransient, TMI, Statistics, InteractiveUtils, Interpolations, PlutoUI, Plots

# ╔═╡ 00d141f0-10b8-458d-b1d6-4caba8193542
md"""
# Vintage Test 
This test demonstrates how to compute a "vintage distribution" - this is the percentage of waters, at every grid cell, sourced from the surface during some time period. In general, this is calculated by solving 

$$\frac{\partial c}{\partial t} = Lc$$

where $\mathbf{c}$ is the concentration of a tracer, $\mathbf{L}$ is the tracer transport operator, and $t$ is time. The $\mathbf{L}$ matrix does not affect the surface condition, so whatever surface condition/initial condition is provided will be sustained throughout the time-stepping. 
"""

# ╔═╡ 88606362-2260-4726-825b-707b030bb09b
begin 
	TMIversion = "modern_90x45x33_GH10_GH12"
	A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion,compute_lu=false);
end

# ╔═╡ 860ff5ea-3fbc-4c96-a453-645f5f6ef5a6
md"""
Read in the `.mat` file with the pre-computed step-response. This is a vector of fields, where each field shows the response of the global interior ocean to constant surface conditions at increasing timesteps. 

Move the sliders to move through time and depth, and see which regions of the ocean have surface waters propagated them to faster.
"""

# ╔═╡ 16fad661-5a84-463f-b5e1-0201a63d43c4
begin
	Δ, τ = read_stepresponse()
	lat = Δ[1].γ.lat
	lon = Δ[1].γ.lon
	depth = Δ[1].γ.depth
end

# ╔═╡ 3b84d7e8-6d71-4c47-887a-b27c49d73f7f
@bind di Slider(depth, show_value = true) 

# ╔═╡ 0f3a9f78-5c9b-4e70-9c1f-94396fe31ea3
@bind τi Slider(τ, show_value = true) 

# ╔═╡ a2d9f5d7-e10a-4461-a1f9-62d6c55e5021
begin
	depth_index = findall(x->x == di, depth)[1]
	τ_index = findall(x-> x== τi, τ)[1]
	contourf(lon, lat, Δ[τ_index].tracer[:, :, depth_index]', title = "Step response @ Depth = " * string(di) *", Time = " * string(τi), colorbar_title = "% surface-sourced waters")
end

# ╔═╡ 39debc60-7253-43e0-9db7-f6673909a0ac
md"""
Now we attempt to compute the % of waters in the modern ocean from a specified pre-modern time period. This can be computed through the `vintagedistribution` function, which either reads in a precomputed `.mat` file or calculates it in Julia (the latter will take longer). 

This computation, for the 4° $\times$ 4° TMI, is as follows

Dye the surface of the ocean by allocating a vector $\mathbf{b}$ of length 2806 (of all surface points, which are wet/ocean?). Turn this into a boundary condition by multiplying it by $\mathbf{B}$, a matrix that maps a surface boundary condition to a Dirichlet boundary condition. $\mathbf{B}$ is 74064 $\times$ 2806, where 74064 is the number of wet interior points. We will call the resulting initial condition $\mathbf{c_0}$

$$\mathbf{c_0} = \mathbf{B} \mathbf{b}$$

Then, solve the differential equation 

$$\frac{\partial \mathbf{c}}{\partial t} = \mathbf{L ⋅ c}(t)$$ 

for all relevant times $t$. $\mathbf{L}$ is the discrete time version of the tracer transport operator $\mathcal{L}$. We are provided the initial condition $\mathbf{c_0}$. The matrix $\mathbf{L}$ has 0s where the initial condition is, 

If we want to compute the percentage of waters in modern ocean from 2015-2020, we simply compute the percentage of waters sourced from the surface between 2015-modern, and the percentage of waters sourced from the surface between 2020-modern. The difference of these two will provide the percentage of waters sourced between 2015-2020. 
"""

# ╔═╡ 33ed7f07-0488-48ca-ac25-87370999672b
@time g = vintagedistribution(TMIversion,γ,L,B,2015,2020)

# ╔═╡ 53aa04e3-352d-4669-b921-ffb53078db27
md"""
For random interior locations, we can compute the percentage of modern-day waters sourced from the surface between 2015-2020 as follows. 
"""

# ╔═╡ 57893e28-7d78-4b7a-9a01-5c9e8eb75371
begin
	N = 2
	# get random locations that are wet (ocean)
	locs = Vector{Tuple{Float64,Float64,Float64}}(undef,N)
	[locs[i] = wetlocation(γ) for i in eachindex(locs)]
	
	# get weighted interpolation indices
	wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
	[wis[i] = interpindex(locs[i],γ) for i in 1:N]
	
	y1 = observe(g,wis,γ)
	#y2 = observe(g2,wis,γ)
	println("Location 1 = " * string(locs[1]) * ", " * string(round(y1[1] * 100, sigdigits = 3)) * "% waters from 2015-2020")
	println("Location 1 = " * string(locs[2]) * ", " * string(round(y1[2] * 100, sigdigits = 3)) * "% waters from 2015-2020")
end

# ╔═╡ 503b6c33-63f3-405d-8b24-dcf37878e464
@bind di2 Slider(depth, show_value = true) 

# ╔═╡ bfc41258-8422-482a-9042-2d1f9e3f866c
begin
	depth_index2 = findall(x->x == di2, depth)[1]
	contourf(lon, lat, g.tracer[:, :, depth_index2]', title = "% of waters in modern ocean from 2015-2020 @ Depth = " * string(di), colorbar_title = "%")
end

# ╔═╡ Cell order:
# ╟─00d141f0-10b8-458d-b1d6-4caba8193542
# ╠═50b8ca96-1c49-4690-a0ae-5f8cacaea429
# ╠═a088c589-906d-4d18-b9d9-44dab6c7fbcd
# ╠═88606362-2260-4726-825b-707b030bb09b
# ╟─860ff5ea-3fbc-4c96-a453-645f5f6ef5a6
# ╠═16fad661-5a84-463f-b5e1-0201a63d43c4
# ╠═3b84d7e8-6d71-4c47-887a-b27c49d73f7f
# ╠═0f3a9f78-5c9b-4e70-9c1f-94396fe31ea3
# ╠═a2d9f5d7-e10a-4461-a1f9-62d6c55e5021
# ╟─39debc60-7253-43e0-9db7-f6673909a0ac
# ╠═33ed7f07-0488-48ca-ac25-87370999672b
# ╟─53aa04e3-352d-4669-b921-ffb53078db27
# ╠═57893e28-7d78-4b7a-9a01-5c9e8eb75371
# ╠═503b6c33-63f3-405d-8b24-dcf37878e464
# ╠═bfc41258-8422-482a-9042-2d1f9e3f866c
