# Required Inputs

As briefly described in the [Quick Start](@ref) there are three required inputs to the analysis function: a `DuctedRotor` object, comprised of duct and center body coordinates, a `Rotor` object, and a `PanelingConstants` object; an `OperatingPoint` or vector of operating point objects; and a `ReferenceParameters` object.

## DuctedRotor

The `DuctedRotor` input contains all the information related to the geometry of the ducted rotor system and how that geometry will be handled in the preprocessing stage of the analysis.

```@docs
DuctAPE.DuctedRotor
```

### Duct and Center Body Coordinates

The duct and center body coordinates are input into the `DuctedRotor` type in the order of duct then center body.
Both sets of coordinates must be given in a counter-clockwise order.
For the duct, that means the coordinates start at the inside trailing edge, and continue clockwise around the leading edge and end at the outside trailing edge.
For the center body, that means the coordintes start at the leading edge and end at the trailing edge.
In addition, the coordinates are expected to be provided in columns of the format [z r], where z are the axial coordinates and r are the radial coordinates (analogous to x and y common for 2D airfoils).
There are several checks in place when the `DuctedRotor` object is constructed to attempt to ensure the formats are correct, but it would be best to input things correctly in the first place.

### Rotor

The `Rotor` object contains mostly information about the rotor blade geometry, but also includes details used for the blade section lift and drag polar lookups.
For the most part, the inputs are pretty straightforward and the docstring below is sufficient to understand what they are.
In general, the blade geometry should be defined from hub to tip, therefore the radial locations should be monotonically increasing.
Care should be taken to note which inputs are dimensional and which are non-dimensional.
It should also be noted that the rotor inputs are _always_ interpolated as part of the analysis preprocess.

```@docs
DuctAPE.Rotor
```

#### Airfoils

The input with the most flexibility in the `Rotor` object is the airfoils input.
These "airfoils" are really the methods for looking up the lift and drag polar values given inputs such as angle of attack, Reynolds number, Mach number, etc.
Airfoil types for DuctAPE are currently contained in the [C$^\textrm{4}$Blade](@ref "C$^\textrm{4}$Blade [[C](#)ascade [C](#)ompatible [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/)]") sub-module of DuctAPE which is exported as `c4b` and also contains the various airfoil evaluation functions used for the blade element lookups.
The available airfoil types include all the airfoil types from CCBlade, as well as `DFDCairfoil` which is an [XROTOR](https://web.mit.edu/drela/Public/web/xrotor/)-like parametric cascade polar used in DFDC.
In addition there are untested cascade types with similar structure to CCBlades airfoil types called `DTCascade`.
Furthermore, there is an experimental actuator disk model implemented via the `ADM` airfoil type in C4Blade.
 Any airfoil object used needs to be able to be deconstructed into a vector for compatiblity with PreallocationTools (see [Precompiled Caches](@ref) below) and then automatedly reconstructed inside the solve.  This process has been generalized for all the C4Blade airfoil types, and should work without much issue for other types as long as the airfoil type fields are arrays of numbers.
In addition, if custom user types are added, they must be of one of the C4Blade abstract airfoil types, otherwise DuctAPE will not know which evalutaion function to use and throw an error.

### PanelingConstants

The `PanelingConstants` object contains all the information required to dictate how the preprocess repanling will take place and almost all of the information required to generate the precompiled caches.

The first two inputs dictate the number panels that should be used in repaneling the region of each body from their respective leading edge to the axial position of the first rotor.
Note that even though we have called this region the inlet, it really controls the number of panels ahead of the first rotor's axial position.
The `num_panels` input is a vector describing the number of panels to be linearly spaced between each discrete location in the system, where discrete locations are rotor axial positions as well as the axial position of the body trailing edges (in the case of the duct, this is the interior trailing edge, even though the exterior trailing edge could be further back in the case where a trailing edge gap is present).

The `dte_minus_cbte` is the duct trailing edge axial position minus the center body trailing edge axial position.
If it is zero, the axial positions of duct and center body trailing edge must be exactly equal.

The `num_wake_sheets` controls how many blade elements are interpolated from the `Rotor` input.  The wake sheets are linearly interpolated along the first rotor from center body to duct, so there are one fewer blade elements than wake sheets.
Note that there is no way to manually control how the wake aligns with rotors aft of the first rotor, so those rotor inputs will be interpolated at the points at which the wake sheets align.

The `wake_length` input defaults to 1.0 and we have found that to be a generally reasonable value.
Note that a shorter or longer wake will not lead to a change in number of panels since the last entry in `num_panels` controls that.

```@docs
DuctAPE.PanelingConstants
```

## OperatingPoint

The `OperatingPoint` inputs contains all the information concering the freestream as well as the rotor rotation rate(s).
Note that the pressure and temperature fields in the struct are currently unused.
For analyses with multiple operating points a vector of these objects should be used, even though many of the fields will be repeated.


```@docs
DuctAPE.OperatingPoint
```

## ReferenceParameters

The `ReferenceParameters` are only used in the postprocessing to compute various coefficients.
Usually, the reference radius is simply the rotor tip radius.
The reference velocity may be selected as the freestream velocity, but often may not be, for example in the case of zero freestream velocity (static conditions), the reference velocity should probably be some non-zero value such that the coefficients don't simply return zeros or infinities.

```@docs
DuctAPE.ReferenceParameters
```

