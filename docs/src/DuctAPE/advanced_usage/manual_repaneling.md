## Circumventing the Automated Geometry Re-paneling

It is not advised to circument the automated geometry re-paneling, but if it must be done, the user needs to provide duct, center_body, and wake nodes conforming to compatible geometry formatting.
The best use case for this is to use previously generated geometry or perhaps geometry exported from DFDC.

The process is not simple, but is possible.
You would have to manually run the dispatches of [`precompute_parameters`](@ref "DuctAPE.precompute_parameters") that take in the the repaneled body nodes and wake grid.
These dispatches exist for this purpose, but there is, by design, no convenience functions at this time to aid the user in easily bypassing the automated repaneling.
