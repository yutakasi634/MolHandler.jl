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
    clip_trajectory(query::Union{Integer, Array{T, 1}, OrdinalRange}, trj::Trajectory;
                    query_key = key::Symbol = :frame)::Trajecotory where T <: Integer

Clip a part of trajectory you specified.
Expected query key is one of the following

- `:frame` : in this case, query specify the frame index.
- `:atom`  : in this case, query specify the atom index.
"""
function clip_trajectory(query::Integer, trj::Trajectory;
                         query_key = query_key::Symbol = :frame)::Trajectory
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

function clip_trajectory(query::Union{Array{T, 1}, OrdinalRange}, trj::Trajectory;
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

"""
    center_of_mass(query::Trajectory;
                   indices = Union{Array, OrdinalRange, Colon} = :,
                   geometric = geometric_flag::Bool = false)::Vector{Coordinate{Float32}}

Calculate the center of mass of trajectory for specified atom indices.
If you set geometric flag is `true`, this function calculate geometric center of mass.
"""
function center_of_mass(trj::Trajectory;
                        indices = indices::Union{Array, OrdinalRange, Colon} = :,
                        geometric = geometric::Bool = false)::Vector{Coordinate{Float32}}
    if !geometric
        mass_vec = map(attributes -> attributes.mass, trj.attributes)
        selected_mass_vec = mass_vec[indices]
        if Nothing ∈ selected_mass_vec
            error("""
                  center_of_mass: geometric flag is false but attributes of trajectory have Nothing value.
                  you can not use mass to calculation.
                  """)
        end
        total_mass = sum(selected_mass_vec)
        weited_sum = reshape(selected_mass_vec, (1, length(selected_mass_vec))) * trj.coordinates[indices, :]
        vec( weited_sum / total_mass)
    else
        selected_coord = trj.coordinates[indices, :]
        vec(sum(selected_coord, dims = 1) / size(selected_coord, 1))
    end
end

function pair_length_matrix(first_coords::Array{Coordinate{T}, 1},
                            second_coords::Array{Coordinate{T}, 1})::Matrix{Coordinate{T}} where T <: Real
    distance.(first_coords, reshape(second_coords, (1, length(second_coords))))
end

function pair_length_matrix(trj::Trajectory;
                            frame_indices = frame_ids::Union{Array, OrdinalRnage, Colon} = :,
                            first_atoms_indices = first_atom_ids::Union{Array, OrdinalRange, Colon} = :,
                            second_atoms_indices = second_atom_ids::Union{Array, OrdinalRnage, Colon} = first_indices)::Vector{Matrix{Coordinate}}
    zip_iterate4frame = zip(each_col(trj.coordinates[first_atom_ids, frame_ids])
                            each_col(trj.coordinates[second_atom_ids, frame_ids]))
    map(coords_vec_pair -> pair_length_matrix(coords_vec_pair...), zip_iterate4frame)
end
