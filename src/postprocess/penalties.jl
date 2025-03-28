"""
    separation_penalty(s_sep, steps, separation_allowance, separation_penalty)

Separation penalty based on boundary layer separation.
Returns a penalty value based on where separation occured relative to the allowed separation region with maximum of the input separation penalty (linearly distributed from zero at begining of region of concern).
"""
function separation_penalty(
    s_sep, steps, separation_allowance, separation_penalty; hardness=1e6
)
    return FLOWMath.ksmax(
        [
            0.0
            FLOWMath.linear(
                [
                    0.0
                    steps[end - separation_allowance]
                ],
                [separation_penalty; 0.0],
                s_sep,
            )
        ],
        hardness,
    )
end
