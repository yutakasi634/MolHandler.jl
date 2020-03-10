"""
    get_frame(frame_idx::Int64, trajectory::Trajectory)
    ::Frame

Return Frame object which correspond to `framd_idx` frame from `trajectory`.
"""
function get_frame(frame_idx::Int64, trj::Trajectory)::Frame
    Frame(trj.coordinates[:,frame_idx], trj.attributes)
end

"""
    get_atom(query::Int64, trajectory::Trajectory)
    ::Vector{Atom}

Return Vector of `query`th Atom object from `trajectory`.
"""
function get_atom(query::Int64, trj::Trajectory)::Vector{Atom}
    map(xyz -> Atom(xyz, trj.attributes[query]), trj.coordinates[query, :])
end

"""
    get_atom(query::Int64, frame::Frame)
    ::Atom

Return `query`th Atom object from `frame`.
"""
function get_atom(query::Int64, frame::Frame)::Atom
    Atom(frame.coordinates[query], frame.attributes[query])
end


"""
    clip_trajectory(query::Union{Integer, Array{T, 1}, OrdinalRange}, trj::Trajectory;
                    query_key::Symbol = :frame)
    ::Trajecotory where T <: Integer

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
                   indices::Union{Array, OrdinalRange, Colon} = :,
                   geometric::Bool = false)
    ::Vector{Coordinate{Float32}}

Calculate the center of mass of trajectory for specified atom indices.
If you set `geometric` is `true`, this function calculate geometric center of mass.
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

"""
    pair_length_matrix(corrds1::Vector{Coordinate}, coords2::Vector{Coordinate})
    ::Matrix{Coordinate}

Calculate distance matrix for all combination between `coords1` and `coord2`.
The row of returned matrix coresspond to `coords1`, and the column correspond to `coords2`.
"""
function pair_length_matrix(first_coords::AbstractArray{Coordinate{T}, 1},
                            second_coords::AbstractArray{Coordinate{T}, 1})::Matrix{T} where T <: Real
    distance.(first_coords, reshape(second_coords, (1, length(second_coords))))
end

"""
    pair_length_matrix(trj::Trajectory;
                       frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                       first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                       second_atom_indices::Union{Array, OrdinalRange, Colon} = atom_indices1)
    ::Vector{Matrix{Coordinate}}

Calculate distance matrix for all combination between `first_atom_indices` and `second_atom_indices`.
This cauculation apply to each frame of trajectory and the result matrices are stored in Vector.
The target frame can be restricted by pass indeces vector or range to `frame_indices`.
"""
function pair_length_matrix(trj::Trajectory;
                            frame_indices::Union{Array, OrdinalRange, Colon} = :,
                            first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                            second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Vector{Matrix{Real}}
    zip_iterate4frame = zip(eachcol(trj.coordinates[first_atom_indices, frame_indices]),
                            eachcol(trj.coordinates[second_atom_indices, frame_indices]))
    map(coords_vec_pair -> pair_length_matrix(coords_vec_pair...), zip_iterate4frame)
end

"""
    contact_bool_matrix(threshold::Real, trj::Trajectory;
                        frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                        first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                        second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Vector{Matrix{Bool}}

Judge contact is formed or not. If the distance between two coordinate is shorter than threshold, contact is considered to be formed. In returned vector of matrices, each matrix correspond to contact matrix of each frame.
You can specify the target frames or atoms by `frame_indices`, `first_atom_indices` or `second_atom_indices`. When you specify the target atoms, the row of matrices corresponds to first_atom_indices and column of matrices corresponds to second_atom_indices.
"""
function contact_bool_matrix(threshold::Real, trj::Trajectory;
                            frame_indices::Union{Array, OrdinalRange, Colon} = :,
                            first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                            second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Vector{Matrix{Bool}}
    length_mat_arr = pair_length_matrix(trj, frame_indices = frame_indices,
                                        first_atom_indices = first_atom_indices,
                                        second_atom_indices = second_atom_indices)
    map(length_mat_arr) do length_matrix
        map(length -> length < threshold, length_matrix)
    end
end

"""
    contact_probability_matrix(threshold::Real, trj::Trajectory;
                               frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                               first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                               second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Matrix{Real}

Calculate contact formation probability over the trajectory. If the distance between two coordinate is shorter than threshold, contact is considered to be formed.
You can specify the target frames or atoms by `frame_indices`, `first_atom_indices` or `second_atom_indices`. When you specify the target atoms, the row of matrices corresponds to first_atom_indices and column of matrices corresponds to second_atom_indices.
"""
function contact_probability_matrix(threshold::Real, trj::Trajectory;
                                    frame_indices::Union{Array, OrdinalRange, Colon} = :,
                                    first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                                    second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Matrix{Real}
    contact_mat_arr = contact_bool_matrix(threshold, trj,
                                          frame_indices = frame_indices,
                                          first_atom_indices = first_atom_indices,
                                          second_atom_indices = second_atom_indices)
    sum(contact_mat_arr) / length(contact_mat_arr)
end
