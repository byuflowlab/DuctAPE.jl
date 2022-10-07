# API Reference

## Public

### Composite Types

There are several native DuctTAPE types that are required for the user to use as inputs.

The first is the freestream object, which defines the various freestream values.

```@docs
DuctTAPE.Freestream
```

Second is the rotor geometry object that contains all the required information do define a rotor as well as allow DuctTAPE to create the CCBlade rotor, section, and operation point inputs.

```@docs
DuctTAPE.RotorGeometry
```

### Functions

The current main function the user will employ from DuctTAPE is `ff2ccb` which is the FLOWFoil to CCBlade coupling function. This function calls several private functions listed below.

```@docs
DuctTAPE.ff2ccb
```

---

## Private

### Auxiliary Coupling Functions

These functions are called by the `ff2ccb` function as part of converting the FLOWFoil outputs into CCBlade inputs.

```@docs
DuctTAPE.extract_ff_geom
DuctTAPE.ff2ccb_velocity
DuctTAPE.generate_ccb_sections
```


### Utility Functions

The following functions are used as various utility functions throughout the code.

```@autodocs
Modules = [DuctTAPE]
Pages   = ["utils.jl"]
```
