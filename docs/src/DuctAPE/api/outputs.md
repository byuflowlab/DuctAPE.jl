# Outputs

The output tuple contains many items.
The [`post_process`](@ref "DuctAPE.post_process") function docstring lists them.
The purpose of showing this function here is not for you to manually run the fuction, but simply rather for you to see what the available outputs are.

```@docs
DuctAPE.post_process
```

## Returning the Pre-process Objects

Sometimes, it may be desireable to return the pre-process objects, including:

- `panels` which is a named tuple containing the body, rotor, and wake panel objects
- `ivb` which are the unit induced velocities on the body panels
- `solve_parameter_tuple` which contains all of the solver parameters
- `blade_elements` which contains all of the blade element geometry and airfoil information
- `linsys` which contains all the linear system objects for the panel method
- `idmaps` which contains all the index mapping used throughout the solve and post-process.

In this case, we can use the `return_inputs` keyword argument when calling the `analyze` function to return a named tuple containing those pre-process objects.

```julia
outs, ins, success_flag = dt.analyze(
    ducted_rotor, operating_point, reference_parameters; return_inputs=true
)
```
