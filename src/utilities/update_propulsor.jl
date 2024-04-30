"""
    update_operating_point!(op_old, op_new)

Overwrites all the values of an OperatingPoint object with another OperatingPoint object's values (or NamedTuple with the same field names).

# Arguments
- `op_old::OperatingPoint` : the OperatingPoint to be overwritten (can also be a NamedTuple with the same field names as an OperatingPoint).
- `op_new::OperatingPoint` : the OperatingPoint values to be used (can also be a NamedTuple with the same field names as an OperatingPoint).
"""
function update_operating_point!(op_old, op_new)
    op_old.Vinf .= op_new.Vinf
    op_old.rhoinf .= op_new.rhoinf
    op_old.muinf .= op_new.muinf
    op_old.asound .= op_new.asound
    op_old.Omega .= op_new.Omega
    return op_old
end
