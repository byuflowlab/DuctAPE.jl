## Multi-Point Analyses

In the case that one wants to run the same geometry at several different operating points, for example: for a range of advance ratios, there is another dispatch of the `analyze` function that takes in an input, `multipoint`, that is a vector of operating points.

```@docs; canonical=false
DuctAPE.analyze(multipoint::AbstractVector{TO},propulsor::Propulsor,options::Options) where TO<:OperatingPoint
```

If we were to continue the [tutorial](@ref "Run Analysis") and run a multi-point analysis on the example geometry given there, it might look something like this:

```julia
# - Advance Ratio Range - #
Js = range(0.0, 2.0; step=0.1)

# - Calculate Vinfs - #
D = 2.0 * rotorstator_parameters.Rtip[1] # rotor diameter
n = RPM / 60.0 # rotation rate in revolutions per second
Vinfs = Js * n * D

# - Set Operating Points - #
ops = [deepcopy(operating_point) for i in 1:length(Vinfs)]
for (iv, v) in enumerate(Vinfs)
    ops[iv].Vinf[] = v
end

# - Run Multi-point Analysis - #
outs_vec, success_flags = dt.analyze(ops, propulsor, dt.set_options(ops))
```

There are a few things to note here.
First, we want to make sure that the operating point objects we put into the input vector are unique instances.
Second, we need to use the dispatch of `set_options` that takes in the operating point vector to set up the right number of things in the background (like convergence flags for each operating point).
Third, the outputs of the analysis are vectors of the same outputs for a single analysis.
