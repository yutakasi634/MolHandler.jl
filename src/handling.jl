"""
    get_frame(frame_idx::Int64, trajectory::Trajectory)::Frame

Return Frame object which correspond to `framd_idx` frame from `trajectory`.
"""
function get_frame(frame_idx::Int64, trj::Trajectory)::Frame
    Frame(trj.coordinates[:,frame_idx], trj.attributes)
end

"""
    get_atom(query::Int64, trajectory::Trajectory)::Vector{Atom}

Return Vector of `query`th Atom object from `trajectory`.
"""
function get_atom(query::Int64, trj::Trajectory)::Vector{Atom}
    map(xyz -> Atom(xyz, trj.attributes[query]), trj.coordinates[query, :])
end

"""
    get_atom(query::Int64, frame::Frame)

Return `query`th Atom object from `frame`.
"""
function get_atom(query::Int64, frame::Frame)::Atom
    Atom(frame.coordinates[query], frame.attributes[query])
end
