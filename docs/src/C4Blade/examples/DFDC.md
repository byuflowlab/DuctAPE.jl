# DFDC Airfoil Type

The DFDC Airfoil type is very similar to the [XROTOR](https://web.mit.edu/drela/Public/web/xrotor/) airfoil type, but includes additions for cascade corrections based on stagger and solidity.
The cascade corrections aren't particularly accurate, but they do apply ballpark effects resulting from high solidity blade sections.
The main benefit to this airfoil type is its simplicity and that the post-stall behavior is already in a format allowing more robust convergence of the DuctAPE solvers.

```@docs
DuctAPE.C4Blade.DFDCairfoil
```

