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
    get_atom(query::Int64, frame::Frame)::Atom

Return `query`th Atom object from `frame`.
"""
function get_atom(query::Int64, frame::Frame)::Atom
    Atom(frame.coordinates[query], frame.attributes[query])
end


"""
    clip_trajectory(query::Union{Integer, Array{T, 1}, OrdinalRange}, trj::Trajectory,
                    query_key::Symbol = :frame)::Trajecotory where T <: Integer

Clip a part of trajectory you specified.
Expected query key is one of the following

    - `:frame` : in this case, query specify the frame index.
    - `:atom`  : in this case, query specify the atom index.
"""
function clip_trajectory(query::Integer, trj::Trajectory,
                         query_key::Symbol = :frame)::Trajectory
    if query_key == :frame
        Trajectory(reshape(trj.coordinates[:, query], (trj.natom, 1)), trj.attributes)
    elseif query_key == :atom
        Trajectory(reshape(trj.coordinates[query, :], (1, trj.nframe)), [trj.attributes[query]])
    else
        error("""
              clip_trajectory: invalid query key :$(query_key) here
              expected query key is one of the following.
              - :frame
              - :atom
              """)
    end
end

    function clip_trajectory(query::Union{Array{T, 1}, OrdinalRange}, trj::Trajectory,
                             query_key::Symbol = :frame)::Trajectory where T <: Integer
    if query_key == :frame
        Trajectory(trj.coordinates[:, query], trj.attributes)
    elseif query_key == :atom
        Trajectory(trj.coordinates[query, :], trj.attributes[query])
    else
        error("""
              clip_trajectory: invalid query key :$(query_key) here
              expected query key is one of the following.
              - :frame
              - :atom
              """)
    end
end

# function center_of_mass(index::Union{Integer, Array, StepRange}, trj::Trajectory)::Vector{Coordinate{Float32}}
# end
