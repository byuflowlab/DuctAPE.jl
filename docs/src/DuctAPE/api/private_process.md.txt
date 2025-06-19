```@contents
Pages = ["private_process.md"]
Depth = 5
```

## Analysis
```@docs
DuctAPE.analyze_multipoint
```

## Process
```@docs
DuctAPE.process
DuctAPE.solve
```

### Internal Solvers
#### ModCSOR
```@docs
DuctAPE.mod_COR_solver
DuctAPE.relax_Gamr_mod!
DuctAPE.relax_gamw_mod!
DuctAPE.update_states!
```

#### CSOR
```@docs
DuctAPE.update_CSOR_residual_values!
DuctAPE.check_CSOR_convergence!
DuctAPE.relax_Gamr!
DuctAPE.relax_gamw!
```

### Residuals

#### ModCSOR
```@docs
DuctAPE.mod_CSOR_residual!
DuctAPE.estimate_CSOR_states!
```

#### CSOR
```@docs
DuctAPE.CSOR_residual!
DuctAPE.compute_CSOR_residual!
```

#### External Solvers
```@docs
DuctAPE.system_residual
DuctAPE.system_residual!
DuctAPE.update_system_residual!
DuctAPE.estimate_states!
```

## Solve Utilities
```@docs
DuctAPE.extract_initial_guess
DuctAPE.extract_state_variables
```


