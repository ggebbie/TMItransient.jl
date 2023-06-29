# What is this? 
Notebooks emulating what happens in `runtests.jl` with some deeper explanations 

# How I Made this 
Assume whoever is using it has a Julia installation with some version of Pluto. We'll have a special version of the environment/manifest that will be activated for the notebooks. 

to handle the Manifest/Project files, I usually work just in terminal 
- in terminal start up julia, and create a new environment in the `notebooks` folder 
- add `TMI` through git
- add `TMItransient` by setting the working directory to where TMItransient is on your machine, then just add it 
- add anything else that you want (Plots, InteractiveUtils) *that is specific to the notebook examples*
  - note: if I just install Plots, it won't precompile (missing WebSocket?), but once I add PlotlyJS, it does install 
- should Pluto be added to this Project/Manifest? Maybe...

to create a new notebook
- in terminal 
``` 
julia
import Pluto; Pluto.run()
```
- *name the new notebook*, make sure you name it before trying to access the packages, otherwise it's in the `.julia` folder and has a `.toml` files in `tmp` 
- in the notebook, add the line `import Pkg; Pkg.activate(".")` so that we activate the environem

