# Analysis

Since the remaining API references here involve the inputs and outputs of the analysis function, we begin by looking at the `analyze` function.
In general, there are two main dispatches that will be used in the vast majority of use cases.
The first is for a single operating point, and the second is for multiple operating points.
They are identical, except for the multiple operating point dispatch requiring a vector of operating points and returning a vector of output objects.

```@docs
DuctAPE.analyze(ducted_rotor::DuctedRotor,operating_point::OperatingPoint,reference_parameters::ReferenceParameters)
DuctAPE.analyze(ducted_rotor::DuctedRotor,operating_point::AbstractVector{OperatingPoint},reference_parameters::ReferenceParameters)
```


