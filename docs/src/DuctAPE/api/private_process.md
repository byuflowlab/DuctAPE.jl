## Analysis
```@docs
DuctAPE.analyze_multipoint
```

### Process
```@docs
DuctAPE.process
DuctAPE.solve
```

#### Residuals

##### CSOR
```@docs
DuctAPE.CSOR_residual!
DuctAPE.compute_CSOR_residual!
DuctAPE.relax_Gamr!
DuctAPE.relax_gamw!
DuctAPE.apply_relaxation_schedule
DuctAPE.update_CSOR_residual_values!
DuctAPE.check_CSOR_convergence!
```

##### External Solvers
```@docs
DuctAPE.system_residual
DuctAPE.system_residual!
DuctAPE.update_system_residual!
DuctAPE.estimate_states!
```

#### Solve Utilities
```@docs
DuctAPE.extract_initial_guess
DuctAPE.extract_state_variables
```


