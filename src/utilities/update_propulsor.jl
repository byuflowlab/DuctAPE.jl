"""
"""
function update_operating_point!(op_old, op_new)
    op_old.Vinf .= op_new.Vinf
    op_old.rhoinf .= op_new.rhoinf
    op_old.muinf .= op_new.muinf
    op_old.asound .= op_new.asound
    op_old.Omega .= op_new.Omega
    return op_old
end
