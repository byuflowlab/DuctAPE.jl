## Available Outputs

The output tuple contains many items.
The [`post_process`](@ref "DuctAPE.post_process") function docstring lists them.

### Returning the Pre-process Objects

Sometimes, it may be desireable to return the pre-process objects, including:

- `panels` which is a named tuple containing the body, rotor, and wake panel objects
- `ivb` which are the unit induced velocities on the body panels
- `solve_parameter_tuple` which contains all of the solver parameters
- `blade_elements` which contains all of the blade element geometry and airfoil information
- `linsys` which contains all the linear system objects for the panel method
- `idmaps` which contains all the index mapping used throughout the solve and post-process.

In this case, we can use the `return_inputs` keyword argument when calling the `analyze` function to return a named tuple containing those pre-process objects.

```julia
outs, ins, success_flag = dt.analyze(propulsor; return_inputs=true)
```
