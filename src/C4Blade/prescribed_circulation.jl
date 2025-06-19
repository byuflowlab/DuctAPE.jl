"""
    ADM
    
# Fields
- `prescribed_circulation::Float=0.0` : Prescribed circulation strength
- `prescribed_source_strength::Float=0.0` : Prescribed source panel strength
"""
@kwdef struct ADM{TFG,TFS}
    prescribed_circulation::TFG = 0.0
    prescribed_source_strength::TFS = 0.0
end
