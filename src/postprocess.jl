#=
Post Processing Types and Functions

Authors: Judd Mehr,
=#

"""
"""
struct Forces{TT;TQ}
    thrust::TT
    torque::TQ
end
