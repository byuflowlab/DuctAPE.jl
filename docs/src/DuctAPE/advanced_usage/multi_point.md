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
outs_vec, success_flags = DuctAPE.analyze(ops, propulsor, DuctAPE.set_options(ops))
```

There are a few things to note here.
First, we want to make sure that the operating point objects we put into the input vector are unique instances.
Second, we need to use the dispatch of `set_options` that takes in the operating point vector to set up the right number of things in the background (like convergence flags for each operating point).
Third, the outputs of the analysis are vectors of the same outputs for a single analysis.


## Advanced Options for Multi-point analyses

For using advanced options in multi-point analyses, there are various changes that need to be made to avoid run-time errors.
Here is an example for setting options with the CSOR solver.

```julia
# number of operating points to analyze
nop = 3

options = DuctAPE.set_options(;
    solver_options=DuctAPE.CSORSolverOptions(
        converged=fill(false, (1, nop)), # need a convergence flag for each operating point
        iterations=zeros(Int, (1, nop)), # need a iteration count for each operating point
        Vconv=ones(nop), # in this case, we need a reference velocity for each operating point
    ),
    write_outputs=fill(false, nop), # we need to know which of the operating point outputs to write
    outfile=fill("", nop), # we need to include names, even if they won't be used.
    output_tuple_name=fill("outs", nop), # we need to include names, even if they won't be used.
)
```

If using a poly-algorithm with a multi-point solve, then each of the solvers needs to have the multiple `converged` and `iterations` fields for each operating point, and the overall solve type needs to have a `converged` and `iterations` field for each solver and each operating point.

```julia
options = DuctAPE.set_options(;
    solver_options=DuctAPE.ChainSolverOptions(;
        solvers=[ # vector of solvers to use in poly-algorithm
            DuctAPE.NLsolveOptions(;
                algorithm=:anderson,
                atol=1e-12,
                iteration_limit=200
                converged=fill(false, (1,nop)), # flags for each operating point
                iterations=zeros(Int, (1,nop)), # counters for each operating point
            ),
            DuctAPE.MinpackOptions(;
                atol=1e-12,
                iteration_limit=100,
                converged=fill(false, (1,nop)),
                iterations=zeros(Int, (1,nop)),
        ],
        converged=fill(false, (2,nop)), # flags for each solver and each operating point
        iterations=fill(0, (2,nop)), # counts for each solver and each operating point
    ),
)
```
