import Base.Threads

"""
    get_frame(frame_idx::Int64, trajectory::Trajectory)
    ::Frame

Return Frame object which correspond to `framd_idx` frame from `trajectory`.
"""
function get_frame(frame_idx::IntT, trj::Trajectory{RealT})::Frame where {IntT <: Integer, RealT <: Real}
    Frame(trj.coordinates[:,frame_idx], trj.attributes)
end

"""
    get_atom(query::Int64, trajectory::Trajectory)
    ::Vector{Atom}

Return Vector of `query`th Atom object from `trajectory`.
"""
function get_atom(query::IntT, trj::Trajectory)::Vector{Atom} where IntT <: Integer
    map(xyz -> Atom(xyz, trj.attributes[query]), trj.coordinates[query, :])
end

"""
    get_atom(query::Int64, frame::Frame)
    ::Atom

Return `query`th Atom object from `frame`.
"""
function get_atom(query::IntT, frame::Frame)::Atom where IntT <: Integer
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
function clip_trajectory(query::IntT where IntT <: Integer, trj::Trajectory{RealT};
                         query_key = query_key::Symbol = :frame)::Trajectory{RealT} where RealT <: Real
    if query_key == :frame
        Trajectory(reshape(trj.coordinates[:, query], (trj.natom, 1)), trj.attributes)
    elseif query_key == :atom
        Trajectory(reshape(trj.coordinates[query, :], (1, trj.nframe)), [trj.attributes[query]])
    else
        throw(ArgumentError("""
                            clip_trajectory: invalid query key :$(query_key) here
                            expected query key is one of the following.
                            - :frame
                            - :atom
                            """))
    end
end

function clip_trajectory(query::Union{Array{IntT, 1}, OrdinalRange}, trj::Trajectory{RealT};
                         query_key::Symbol = :frame)::Trajectory{RealT} where {IntT <: Integer, RealT <: Real}
    if query_key == :frame
        Trajectory(trj.coordinates[:, query], trj.attributes)
    elseif query_key == :atom
        Trajectory(trj.coordinates[query, :], trj.attributes[query])
    else
        throw(ArgumentError("""
                            clip_trajectory: invalid query key :$(query_key) here
                            expected query key is one of the following.
                            - :frame
                            - :atom
                            """))
    end
end

"""
    center_of_mass(query::Trajectory;
                   frame_indices::Union{Array, OrdinalRange, Colon} = :,
                   atom_indices::Union{Array, OrdinalRange, Colon} = :,
                   geometric::Bool = false)
    ::Vector{Coordinate{Real}}

Calculate the center of mass of trajectory for specified frame and atom indices.
If you set `geometric` is `true`, this function calculate geometric center of mass.
"""
function center_of_mass(trj::Trajectory{RealT};
                        frame_indices::Union{Array, OrdinalRange, Colon} = :,
                        atom_indices::Union{Array, OrdinalRange, Colon} = :,
                        geometric::Bool = false)::Vector{Coordinate{RealT}} where RealT <: Real
    if !geometric
        mass_vec = map(attributes -> attributes.mass, trj.attributes)
        selected_mass_vec = view(mass_vec, atom_indices)
        if nothing ∈ selected_mass_vec || length(selected_mass_vec) == 0
            throw(ArgumentError("""
                                center_of_mass: geometric flag is false but attributes of trajectory have Nothing value.
                                you can not use mass to calculation.
                                """))
        end
        total_mass = sum(selected_mass_vec)
        weited_sum = reshape(selected_mass_vec, (1, length(selected_mass_vec))) * view(trj.coordinates, atom_indices, frame_indices)
        vec(weited_sum / total_mass)
    else
        selected_coord = view(trj.coordinates, atom_indices, frame_indices)
        vec(sum(selected_coord, dims = 1) / size(selected_coord, 1))
    end
end

"""
    pair_length_matrix(corrds1::Vector{Coordinate}, coords2::Vector{Coordinate})
    ::Matrix{Coordinate}

Calculate distance matrix for all combination between `coords1` and `coord2`.
The row of returned matrix coresspond to `coords1`, and the column correspond to `coords2`.
"""
function pair_length_matrix(first_coords::AbstractArray{Coordinate{RealT}, 1},
                            second_coords::AbstractArray{Coordinate{RealT}, 1})::Matrix{RealT} where RealT <: Real
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
function pair_length_matrix(trj::Trajectory{RealT};
                            frame_indices::Union{Array, OrdinalRange, Colon} = :,
                            first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                            second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Vector{Matrix{RealT}} where RealT <: Real
    zip_iterate4frame = zip(eachcol(trj.coordinates[first_atom_indices, frame_indices]),
                            eachcol(trj.coordinates[second_atom_indices, frame_indices]))
    map(coords_vec_pair -> pair_length_matrix(coords_vec_pair...), zip_iterate4frame)
end

"""
    pair_length_matrix_parallel(trj::Trajectory;
                                frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                                first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                                second_atom_indices::Union{Array, OrdinalRange, Colon} = atom_indices1)
    ::Vector{Matrix{Coordinate}}

Multi-threads version of pair_length_matrix. If you set available threads number to `Threads.nthreads()`, this function would faster than non-parallel version.
"""
function pair_length_matrix_parallel(trj::Trajectory{RealT};
                                     frame_indices::Union{Array, OrdinalRange, Colon} = :,
                                     first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                                     second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Vector{Matrix{RealT}} where RealT <: Real
    first_target_coords  = trj.coordinates[first_atom_indices,  frame_indices]
    second_target_coords = trj.coordinates[second_atom_indices, frame_indices]
    target_frame_num = size(first_target_coords, 2)
    return_vec = Vector{Matrix{Float32}}(undef, target_frame_num)
    Threads.@threads for frame_idx in 1:target_frame_num
        return_vec[frame_idx] =
            pair_length_matrix(first_target_coords[:, frame_idx], second_target_coords[:, frame_idx])
    end
    return_vec
end

"""
    contact_bool_matrix(threshold::Real, trj::Trajectory;
                        frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                        first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                        second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Vector{Matrix{Bool}}

Judge contact is formed or not. If the distance between two coordinate is shorter than threshold, contact is considered to be formed. In returned matrix, the row of matrices corresponds to first_atom_indices and column of matrices corresponds to second_atom_indices.
"""
function contact_bool_matrix(threshold::RealT,
                             first_coords::AbstractArray{Coordinate{RealT}, 1},
                             second_coords::AbstractArray{Coordinate{RealT}, 1})::Matrix{Bool} where RealT <: Real
    length_matrix = pair_length_matrix(first_coords, second_coords)
    map(length -> length < threshold, length_matrix)
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
function contact_bool_matrix(threshold::RealT1, trj::Trajectory{RealT2};
                            frame_indices::Union{Array, OrdinalRange, Colon} = :,
                            first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                            second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Vector{Matrix{Bool}} where {RealT1 <: Real, RealT2 <: Real}
    length_mat_arr = pair_length_matrix(trj, frame_indices = frame_indices,
                                        first_atom_indices = first_atom_indices,
                                        second_atom_indices = second_atom_indices)
    map(length_mat_arr) do length_matrix
        map(length -> length < threshold, length_matrix)
    end
end

"""
    contact_bool_matrix_parallel(threshold::Real, trj::Trajectory;
                                 frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                                 first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                                 second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Vector{Matrix{Bool}}

Multi-threads version of contact_bool_matrix. If you set available threads number to `Threads.nthreads()`, this function would faster than non-parallel version.
"""
function contact_bool_matrix_parallel(threshold::RealT1, trj::Trajectory{RealT2};
                                      frame_indices::Union{Array, OrdinalRange, Colon} = :,
                                      first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                                      second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Vector{Matrix{Bool}} where {RealT1 <: Real, RealT2 <: Real}
    length_mat_arr = pair_length_matrix_parallel(trj, frame_indices = frame_indices,
                                                 first_atom_indices = first_atom_indices,
                                                 second_atom_indices = second_atom_indices)
    target_frame_num = length(length_mat_arr)
    return_vec = Vector{Matrix{Float32}}(undef, target_frame_num)
    threads_num = Threads.nthreads()
    iteration_in_thread = Int32(floor(target_frame_num / threads_num))
    Threads.@threads for frame_idx in 1:target_frame_num
        return_vec[frame_idx] = map(length -> length < threshold, length_mat_arr[frame_idx])
    end
    return_vec
end

"""
    contact_count_matrix(threshold::Real, trj::Trajectory;
                         frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                         first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                         second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Matrix{Real}

Calculate contact formation count over the trajectory. If the distance between two coordinate is shorter than threshold, contact is considered to be formed.
You can specify the target frames or atoms by `frame_indices`, `first_atom_indices` or `second_atom_indices`. When you specify the target atoms, the row of matrices corresponds to first_atom_indices and column of matrices corresponds to second_atom_indices.
"""
function contact_count_matrix(threshold::RealT1, trj::Trajectory{RealT2};
                              frame_indices::Union{Array, OrdinalRange, Colon} = :,
                              first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                              second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Matrix{RealT2} where {RealT1 <: Real, RealT2 <: Real}
    first_target_coords  = trj.coordinates[first_atom_indices,  frame_indices]
    second_target_coords = trj.coordinates[second_atom_indices, frame_indices]
    target_frame_num = size(first_target_coords, 2)
    count_mat = zeros(RealT2, (size(first_target_coords, 1), size(second_target_coords, 1)))
    for frame_idx in 1:target_frame_num
        count_mat +=
            contact_bool_matrix(threshold, first_target_coords[:, frame_idx], second_target_coords[:, frame_idx])
    end
    count_mat
end

"""
    contact_count_matrix_parallel(threshold::Real, trj::Trajectory;
                                  frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                                  first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                                  second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Matrix{Real}

Multi-threads version of contact_count_matrix. If you set available threads number to `Threads.nthreads()`, this function would faster than non-parallel version.
"""
function contact_count_matrix_parallel(threshold::RealT1, trj::Trajectory{RealT2};
                                       frame_indices::Union{Array, OrdinalRange, Colon} = :,
                                       first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                                       second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Matrix{RealT2} where {RealT1 <: Real, RealT2 <: Real}
    first_target_coords  = trj.coordinates[first_atom_indices,  frame_indices]
    second_target_coords = trj.coordinates[second_atom_indices, frame_indices]
    target_frame_num = size(first_target_coords, 2)
    first_atoms_len = size(first_target_coords, 1)
    second_atoms_len = size(second_target_coords, 1)
    count_mat_arr =
        [zeros(RealT2, first_atoms_len, second_atoms_len) for thread in 1:Threads.nthreads()]
    Threads.@threads for frame_idx in 1:target_frame_num
        count_mat_arr[Threads.threadid()] +=
            contact_bool_matrix(threshold,
                                first_target_coords[:, frame_idx], second_target_coords[:, frame_idx])
    end
    sum(count_mat_arr)
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
function contact_probability_matrix(threshold::RealT1, trj::Trajectory{RealT2};
                                    frame_indices::Union{Array, OrdinalRange, Colon} = :,
                                    first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                                    second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Matrix{RealT2} where {RealT1 <: Real, RealT2 <: Real}
    count_mat= contact_count_matrix(threshold, trj,
                                    frame_indices       = frame_indices,
                                    first_atom_indices  = first_atom_indices,
                                    second_atom_indices = second_atom_indices)
    count_mat / size(trj.coordinates[:, frame_indices], 2)
end

"""
    contact_probability_matrix_parallel(threshold::Real, trj::Trajectory;
                                        frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                                        first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                                        second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Matrix{Real}

Multi-threads version of contact_probability_matrix. If you set available threads number to `Threads.nthreads()`, this function would faster than non-parallel version.
"""
function contact_probability_matrix_parallel(threshold::RealT1, trj::Trajectory{RealT2};
                                             frame_indices::Union{Array, OrdinalRange, Colon} = :,
                                             first_atom_indices::Union{Array, OrdinalRange, Colon} = :,
                                             second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)::Matrix{RealT2} where {RealT1 <: Real, RealT2 <: Real}
    count_mat_arr = contact_count_matrix_parallel(threshold, trj,
                                                  frame_indices       = frame_indices,
                                                  first_atom_indices  = first_atom_indices,
                                                  second_atom_indices = second_atom_indices)
    count_mat_arr / size(trj.coordinates[:,frame_indices], 2)
end

"""
    radius_of_gyration(trj::Trajectory;
                       frame_indices::Union{Array, OrdinalRange, Colon} = :,
                       atom_indices ::Union{Array, OrdinalRange, Colon} = :,
                       geometric::Bool = false)
    ::Vector{Real}

Calculate radius of gyration over the trajectory.
If you set `geometric` is `true`, this function calculate radius of gyration without mass info.
"""
function radius_of_gyration(trj::Trajectory{RealT};
                            frame_indices::Union{Array, OrdinalRange, Colon} = :,
                            atom_indices ::Union{Array, OrdinalRange, Colon} = :,
                            geometric::Bool = false)::Vector{RealT} where RealT <: Real
    com = center_of_mass(trj,
                         frame_indices = frame_indices,
                         atom_indices  = atom_indices,
                         geometric     = geometric)

    target_coordinates = view(trj.coordinates, atom_indices, frame_indices)
    target_frame_num = length(com)
    return_vec = Vector{RealT}(undef, target_frame_num)

    dist_from_com_vec = norm.(target_coordinates .- permutedims(com))

    if !geometric
        mass_vec = map(attributes -> attributes.mass, view(trj.attributes, atom_indices))
        mass_sum = sum(mass_vec)
        for frame_idx in 1:target_frame_num
            return_vec[frame_idx] = sqrt(sum(view(dist_from_com_vec, :, frame_idx).^2 .* mass_vec) / mass_sum)
        end
    else
        target_atom_num = size(target_coordinates)[1]
        for frame_idx in 1:target_frame_num
            return_vec[frame_idx] = sqrt(sum(view(dist_from_com_vec, :, frame_idx).^2) / target_atom_num)
        end
    end
    return_vec
end

"""
    atom_mass(name::AbstractString)
    ::Float32

Return mass of the `name` atom.
The kind of each atom is judged by PDB format atom names rule.
https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf
"""
function atom_mass(name::AbstractString)::Float32
    if occursin(r"^(H|[1-4]H)", name)
        atomname2mass["H"]
    elseif occursin(r"^C", name)
        atomname2mass["C"]
    elseif occursin(r"^O", name)
        atomname2mass["O"]
    elseif occursin(r"^N($|[^A]+)", name)
        atomname2mass["N"]
    elseif occursin(r"^S", name)
        atomname2mass["S"]
    elseif occursin(r"^P", name)
        atomname2mass["P"]
    else
        atomname2mass[name]
    end
end

"""
    residue_mass(name::AbstractString)
    ::Float32

Return mass of the `name` residue.
The kind of each residue is judged by amino acid 3-letter abbreviation.
"""
function residue_mass(name::AbstractString)::Float32
    resname2mass[name]
end
