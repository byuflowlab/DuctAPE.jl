# Cascade Types

!!! warning "Warning"
    Cascade types are currently in development and not ready for general use.

Cascade types are defined similarly to CCBlade airfoil types.
Instead of angle of attack, however, cascade types take in both inflow and stagger angles.
In addition, cascade types are dependent on local solidity.
Currently only one such method is implemented, taking in inflow, stagger, solidity, Reynolds number, and Mach number.  Thus the lift and drag fields for this type are 5D arrays.
Note that since the chord and twist are already fixed inputs to DuctAPE, the solidity and stagger angles are already known a priori.  In many cases, the CCBlade-like airfoil types are sufficient not only for airfoils, but also higher solidity cases.

```@docs
DuctAPE.C4Blade.InReStSoMaCAS
```
