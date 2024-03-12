"""
- `mcrit::Float` : critical Mach number
"""
struct ADM{TFG,TFS}
    prescribed_circulation::TFG
    prescribed_source_strength::TFS
end

function ADM(;
    prescribed_circulation=0.0,
    prescribed_source_strength=0.0,
)
    return ADM(prescribed_circulation, prescribed_source_strength)
end
