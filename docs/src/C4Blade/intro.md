# C$^\textrm{4}$Blade [[C](#)ascade [C](#)ompatible [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/)]

C$^4$Blade is a DuctAPE submodule containing a modified version of [CCBlade](https://flow.byu.edu/CCBlade.jl/stable/) for use in DuctAPE as well as several additional helper functions related to blade section polars.


## Airfoil Types

|Type|Inputs|Status|
|---|---|---|
|[AlphaAF](@ref "DuctAPE.C4Blade.AlphaAF")|angle of attack|✅|
|[AlphaReAF](@ref "DuctAPE.C4Blade.AlphaReAF")|angle of attack, Reynolds|✅|
|[AlphaMachAF](@ref "DuctAPE.C4Blade.AlphaMachAF")|angle of attack, Mach|✅|
|[AlphaReMachAF](@ref "DuctAPE.C4Blade.AlphaReMachAF")|angle of attack, Reynolds, Mach|✅|
|[DFDCairfoil](@ref "DuctAPE.C4Blade.DFDCairfoil")|angle of attack, Reynolds, Mach, solidity, stagger|✅|
|[ADM](@ref "DuctAPE.C4Blade.ADM")|none|🚧|
|[InReStSoMaCAS](@ref "DuctAPE.C4Blade.InReStSoMaCAS")|inflow angle, Reynolds, stagger, solidity, Mach|🚧|

Key:
- ✅ Implemented
- 🚧 Under Development


