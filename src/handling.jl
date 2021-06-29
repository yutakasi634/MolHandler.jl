import Base.Threads
import LinearAlgebra

"""
    get_frame(frame_idx::Int64, trajectory::Trajectory)
    ::Frame

Return Frame object which correspond to `framd_idx` frame from `trajectory`.
"""
function get_frame(frame_idx::IntT, trj::Trajectory{RealT}
    )::Frame where {IntT <: Integer, RealT <: Real}

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
    clip_trajectory(query::Union{Integer, Vector{IntT}, OrdinalRange}, trj::Trajectory;
                    query_key::Symbol = :frame)
    ::Trajecotory where IntT <: Integer

Clip a part of trajectory you specified.
Expected query key is one of the following

- `:frame` : in this case, query specify the frame index.
- `:atom`  : in this case, query specify the atom index.
"""
function clip_trajectory(query::IntT, trj::Trajectory{RealT};
    query_key = query_key::Symbol = :frame
    )::Trajectory{RealT} where {IntT <: Integer, RealT <: Real}

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

function clip_trajectory(query::Union{Vector{IntT}, OrdinalRange}, trj::Trajectory{RealT};
    query_key::Symbol = :frame
    )::Trajectory{RealT} where {IntT <: Integer, RealT <: Real}

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
    geometric_center_of_mass(coordinates::Vector{Coordinate};
                             atom_indices::Union{Vector, OrdinalRange, Colon} = :)
    ::Coordinate

Calculate the geometric center of mass of Coordinate vector for specified atom indices.

> For trajectory case, geometric center of mass is calculated by `center_of_mass` method with geometric flag.
"""
function geometric_center_of_mass(coordinates::Vector{<:Coordinate{RealT}};
    atom_indices::Union{Vector, OrdinalRange, Colon} = :)::Coordinate{RealT} where RealT <: Real

    selected_coord = view(coordinates, atom_indices)
    sum(selected_coord) / length(selected_coord)
end

"""
    center_of_mass(query::Trajectory;
                   frame_indices::Union{Vector, OrdinalRange, Colon} = :,
                   atom_indices::Union{Vector, OrdinalRange, Colon} = :,
                   geometric::Bool = false)
    ::Vector{Coordinate{Real}}

Calculate the center of mass of trajectory for specified frame and atom indices.
If you set `geometric` is `true`, this function calculate geometric center of mass.
"""
function center_of_mass(trj::Trajectory{RealT};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    atom_indices::Union{Vector, OrdinalRange, Colon} = :,
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
    ::Matrix{Real}

Calculate distance matrix for all combination between `coords1` and `coord2`.
The row of returned matrix coresspond to `coords1`, and the column correspond to `coords2`.
"""
function pair_length_matrix(first_coords::ArrayT1, second_coords::ArrayT2
    )::Matrix{<:Real} where {ArrayT1 <: AbstractArray{<:Coordinate{<:Real}, 1},
                             ArrayT2 <: AbstractArray{<:Coordinate{<:Real}, 1}}

    distance.(first_coords, reshape(second_coords, (1, length(second_coords))))
end

"""
     pair_length_matrix_pbc(corrds1::Vector{Coordinate}, coords2::Vector{Coordinate},
                            box_size::Coordinate)
    ::Matrix{Real}

Calculate distance matrix for all combination between `coords1` and `coord2` considering periodic boundary condition.
The box size information used for calculation is specified by `box_size` argument.
The row of returned matrix coresspond to `coords1`, and the column correspond to `coords2`.
"""
function pair_length_matrix_pbc(first_coords::ArrayT1, second_coords::ArrayT2,
    box_size::Coordinate{<:Real}
    )::Matrix{<:Real} where {ArrayT1 <: AbstractArray{<:Coordinate{<:Real}, 1},
                             ArrayT2 <: AbstractArray{<:Coordinate{<:Real}, 1}}

    distance_pbc.(first_coords, reshape(second_coords, (1, length(second_coords))), box_size)
end

"""
    pair_length_matrix(trj::Trajectory;
                       frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                       first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                       second_atom_indices::Union{Vector, OrdinalRange, Colon} = atom_indices1)
    ::Vector{Matrix{Coordinate}}

Calculate distance matrix for all combination between `first_atom_indices` and `second_atom_indices`.
This cauculation apply to each frame of trajectory and the result matrices are stored in Vector.
The target frame can be restricted by pass indeces vector or range to `frame_indices`.
"""
function pair_length_matrix(trj::Trajectory{<:Real};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Vector{<:Matrix{<:Real}}

    zip_iterate4frame = zip(eachcol(view(trj.coordinates, first_atom_indices, frame_indices)),
                            eachcol(view(trj.coordinates, second_atom_indices, frame_indices)))
    map(coords_vec_pair -> pair_length_matrix(coords_vec_pair...), zip_iterate4frame)
end

"""
    pair_length_matrix_pbc(trj::Trajectory, box_size::Coordinate;
                           frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                           first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                           second_atom_indices::Union{Vector, OrdinalRange, Colon} = atom_indices1)
    ::Vector{Matrix{Coordinate}}

Calculate distance matrix for all combinations between `first_atom_indices` and `second_atom_indices` considering periodic boundary condition.
The box size information used for calculation is specified by `box_size` argument.
This cauculation apply to each frame of trajectory and the result matrices are stored in Vector.
The target frame can be restricted by pass indeces vector or range to `frame_indices`.
"""
function pair_length_matrix_pbc(trj::Trajectory{<:Real}, box_size::Coordinate{<:Real};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Vector{<:Matrix{<:Real}}

    zip_iterate4frame = zip(eachcol(view(trj.coordinates, first_atom_indices, frame_indices)),
                            eachcol(view(trj.coordinates, second_atom_indices, frame_indices)))
    map(coords_vec_pair -> pair_length_matrix_pbc(coords_vec_pair..., box_size),
        zip_iterate4frame)
end

"""
    pair_length_matrix_parallel(trj::Trajectory;
                                frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                                first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                                second_atom_indices::Union{Vector, OrdinalRange, Colon} = atom_indices1)
    ::Vector{Matrix{Coordinate}}

Multi-threads version of pair_length_matrix. If you set available threads number to `Threads.nthreads()`, this function would faster than non-parallel version.
"""
function pair_length_matrix_parallel(trj::Trajectory{RealT};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Vector{Matrix{RealT}} where RealT <: Real

    first_target_coords  = view(trj.coordinates, first_atom_indices,  frame_indices)
    second_target_coords = view(trj.coordinates, second_atom_indices, frame_indices)
    target_frame_num = size(first_target_coords, 2)
    return_vec = Vector{Matrix{Float32}}(undef, target_frame_num)
    Threads.@threads for frame_idx in 1:target_frame_num
        return_vec[frame_idx] =
            pair_length_matrix(view(first_target_coords, :, frame_idx), view(second_target_coords, :, frame_idx))
    end
    return_vec
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
function contact_bool_matrix(threshold::RealT,
    first_coords::AbstractArray{<:Coordinate, 1},
    second_coords::AbstractArray{<:Coordinate, 1}
    )::Matrix{Bool} where RealT <: Real

    length_matrix = pair_length_matrix(first_coords, second_coords)
    map(length -> length < threshold, length_matrix)
end

function contact_bool_matrix(threshold::RealT1, trj::Trajectory{RealT2};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Vector{Matrix{Bool}} where {RealT1 <: Real, RealT2 <: Real}

    length_mat_arr = pair_length_matrix(trj, frame_indices = frame_indices,
                                        first_atom_indices = first_atom_indices,
                                        second_atom_indices = second_atom_indices)
    map(length_mat_arr) do length_matrix
        map(length -> length < threshold, length_matrix)
    end
end

"""
    contact_bool_matrix_pbc(threshold::Real, trj::Trajectory, box_size::Coordinate;
                            frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                            first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                            second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Vector{Matrix{Bool}}

Judge contact is formed or not considering periodic boundary condition. If the distance between two coordinate is shorter than threshold, contact is considered to be formed. In returned vector of matrices, each matrix correspond to contact matrix of each frame.
The box size information used for calculation is specified `box_size` argument.
You can specify the target frames or atoms by `frame_indices`, `first_atom_indices` or `second_atom_indices`. When you specify the target atoms, the row of matrices corresponds to first_atom_indices and column of matrices corresponds to second_atom_indices.
"""
function contact_bool_matrix_pbc(threshold::RealT, trj::Trajectory{<:Real},
    box_size::Coordinate{<:Real};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Vector{Matrix{Bool}} where RealT <: Real

    length_mat_arr = pair_length_matrix_pbc(trj, box_size,
                                            frame_indices = frame_indices,
                                            first_atom_indices = first_atom_indices,
                                            second_atom_indices = second_atom_indices)
    map(length_mat_arr) do length_matrix
        map(length -> length < threshold, length_matrix)
    end
end

"""
    contact_bool_matrix_parallel(threshold::Real, trj::Trajectory;
                                 frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                                 first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                                 second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices)
    ::Vector{Matrix{Bool}}

Multi-threads version of contact_bool_matrix. If you set available threads number to `Threads.nthreads()`, this function would faster than non-parallel version.
"""
function contact_bool_matrix_parallel(threshold::RealT1, trj::Trajectory{RealT2};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Vector{Matrix{Bool}} where {RealT1 <: Real, RealT2 <: Real}

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
                         frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                         first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                         second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices)
    ::Matrix{Real}

Calculate contact formation count over the trajectory. If the distance between two coordinate is shorter than threshold, contact is considered to be formed.
You can specify the target frames or atoms by `frame_indices`, `first_atom_indices` or `second_atom_indices`. When you specify the target atoms, the row of matrices corresponds to first_atom_indices and column of matrices corresponds to second_atom_indices.
"""
function contact_count_matrix(threshold::RealT1, trj::Trajectory{RealT2};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Matrix{RealT2} where {RealT1 <: Real, RealT2 <: Real}

    first_target_coords  = view(trj.coordinates, first_atom_indices,  frame_indices)
    second_target_coords = view(trj.coordinates, second_atom_indices, frame_indices)
    target_frame_num = size(first_target_coords, 2)
    count_mat = zeros(RealT2, (size(first_target_coords, 1), size(second_target_coords, 1)))
    for frame_idx in 1:target_frame_num
        count_mat +=
            contact_bool_matrix(threshold, view(first_target_coords, :, frame_idx), view(second_target_coords, :, frame_idx))
    end
    count_mat
end

"""
    contact_count_matrix_parallel(threshold::Real, trj::Trajectory;
                                  frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                                  first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                                  second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices)
    ::Matrix{Real}

Multi-threads version of contact_count_matrix. If you set available threads number to `Threads.nthreads()`, this function would faster than non-parallel version.
"""
function contact_count_matrix_parallel(threshold::RealT1, trj::Trajectory{RealT2};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Matrix{RealT2} where {RealT1 <: Real, RealT2 <: Real}

    first_target_coords  = view(trj.coordinates, first_atom_indices,  frame_indices)
    second_target_coords = view(trj.coordinates, second_atom_indices, frame_indices)
    target_frame_num = size(first_target_coords, 2)
    first_atoms_len = size(first_target_coords, 1)
    second_atoms_len = size(second_target_coords, 1)
    count_mat_arr =
        [zeros(RealT2, first_atoms_len, second_atoms_len) for thread in 1:Threads.nthreads()]
    Threads.@threads for frame_idx in 1:target_frame_num
        count_mat_arr[Threads.threadid()] +=
            contact_bool_matrix(threshold,
                                view(first_target_coords, :, frame_idx), view(second_target_coords, :, frame_idx))
    end
    sum(count_mat_arr)
end

"""
    contact_probability_matrix(threshold::Real, trj::Trajectory;
                               frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                               first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                               second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices)
    ::Matrix{Real}

Calculate contact formation probability over the trajectory. If the distance between two coordinate is shorter than threshold, contact is considered to be formed.
You can specify the target frames or atoms by `frame_indices`, `first_atom_indices` or `second_atom_indices`. When you specify the target atoms, the row of matrices corresponds to first_atom_indices and column of matrices corresponds to second_atom_indices.
"""
function contact_probability_matrix(threshold::RealT1, trj::Trajectory{RealT2};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Matrix{RealT2} where {RealT1 <: Real, RealT2 <: Real}

    count_mat= contact_count_matrix(threshold, trj,
                                    frame_indices       = frame_indices,
                                    first_atom_indices  = first_atom_indices,
                                    second_atom_indices = second_atom_indices)
    count_mat / size(view(trj.coordinates, :, frame_indices), 2)
end

"""
    contact_probability_matrix_parallel(threshold::Real, trj::Trajectory;
                                        frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                                        first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                                        second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices)
    ::Matrix{Real}

Multi-threads version of contact_probability_matrix. If you set available threads number to `Threads.nthreads()`, this function would faster than non-parallel version.
"""
function contact_probability_matrix_parallel(threshold::RealT1, trj::Trajectory{RealT2};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Matrix{RealT2} where {RealT1 <: Real, RealT2 <: Real}

    count_mat_arr = contact_count_matrix_parallel(threshold, trj,
                                                  frame_indices       = frame_indices,
                                                  first_atom_indices  = first_atom_indices,
                                                  second_atom_indices = second_atom_indices)
    count_mat_arr / size(view(trj.coordinates, :,frame_indices), 2)
end

"""
    radius_of_gyration(trj::Trajectory;
                       frame_indices::Union{Vector, OrdinalRange, Colon} = :,
                       atom_indices ::Union{Vector, OrdinalRange, Colon} = :,
                       geometric::Bool = false)
    ::Vector{Real}

Calculate radius of gyration over the trajectory.
If you set `geometric` is `true`, this function calculate radius of gyration without mass info.
"""
function radius_of_gyration(trj::Trajectory{RealT};
             frame_indices::Union{Vector, OrdinalRange, Colon} = :,
             atom_indices ::Union{Vector, OrdinalRange, Colon} = :,
             geometric::Bool = false
             )::Vector{RealT} where RealT <: Real

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
    else target_atom_num = size(target_coordinates)[1]
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
        ATOMNAME2MASS["H"]
    elseif occursin(r"^C($|[^L]+)", name)
        ATOMNAME2MASS["C"]
    elseif occursin(r"^O", name)
        ATOMNAME2MASS["O"]
    elseif occursin(r"^N($|[^A]+)", name)
        ATOMNAME2MASS["N"]
    elseif occursin(r"^S", name)
        ATOMNAME2MASS["S"]
    elseif occursin(r"^P", name)
        ATOMNAME2MASS["P"]
    else
        ATOMNAME2MASS[name]
    end
end

"""
    residue_mass(name::AbstractString)
    ::Float32

Return mass of the `name` residue.
The kind of each residue is judged by amino acid 3-letter abbreviation.
"""
function residue_mass(name::AbstractString)::Float32
    RESNAME2MASS[name]
end

"""
    distance_pdc(first_atom::Coordinate,  second_atom::Coordinate,
                 box_size::Coordinate)
    ::Real

Calculate distance `first_atom` and `second_atom` considering periodic boundary condition.
The box size information used calculatin is specified by `box_size` argument.
"""
function distance_pbc(first_coord::Coordinate{RealT},  second_coord::Coordinate{RealT},
    box_size::Coordinate{<:Real}
    )::RealT where RealT <: Real

    half_box_size = box_size * 0.5
    dist_vec = first_coord  - second_coord
    x = abs(dist_vec.x) < half_box_size.x ? dist_vec.x : box_size.x - abs(dist_vec.x)
    y = abs(dist_vec.y) < half_box_size.y ? dist_vec.y : box_size.y - abs(dist_vec.y)
    z = abs(dist_vec.z) < half_box_size.z ? dist_vec.z : box_size.z - abs(dist_vec.z)
    sqrt(x^2 + y^2 + z^2)
end

# This function get half_box_size, because if process will apply to all trajectory,
# half box size calculation will occur every time frame, this is overhead.
function fix_pbc(coordinates::Vector{<:Coordinate{RealT}}, groupid_vec::Vector{<:Integer},
    box_size::Coordinate{<:Real}, half_box_size::Coordinate{<:Real};
    x::Bool= true, y::Bool = true, z::Bool = true
    )::Vector{<:Coordinate{RealT}} where RealT <: Real

    if length(coordinates) != length(groupid_vec)
        throw(AssertionError("""
                             coordinate vector and groupid vector should be same length.
                             """))
    end

    new_coords = deepcopy(coordinates)
    unique_groupid_vec = unique(groupid_vec)

    for groupid in unique_groupid_vec
        same_group_ids = findall(id->id==groupid, groupid_vec)
        @views sbj_coords = new_coords[same_group_ids[2:end]]
        @views dist2first = sbj_coords .- coordinates[same_group_ids[1]]
        sbj_dist_zip      = zip(sbj_coords, dist2first)

        if x
            for (coord, dist) in sbj_dist_zip
                coord.x = abs(dist.x) < half_box_size.x ? coord.x : coord.x - sign(dist.x) * box_size.x
            end
        end

        if y
            for (coord, dist) in sbj_dist_zip
                coord.y = abs(dist.y) < half_box_size.y ? coord.y : coord.y - sign(dist.y) * box_size.y
            end
        end

        if z
            for (coord, dist) in sbj_dist_zip
                coord.z = abs(dist.z) < half_box_size.z ? coord.z : coord.z - sign(dist.z) * box_size.z
            end
        end
    end
    new_coords
end

"""
    fix_pbc(coordinates::Vector{Coordinate}, groupid_vec::Vector{Integer},
    box_size::Coordinate{<:Real})::Vector{Coordinate}

    Fix the atom group splited by periodic boundary box. This group is defined by groupid_vec which contain each particles group id.
For example, if your system have 3 atom, and atom 1 and 2 are group 1, and atom 3 is group 2, this groupid_vec is [1, 1, 2].
If the group over the boundary of the box, the atoms in the group which separated from first atom of that move to the position where the atom locate without periodic boundary condition.
"""
function fix_pbc(coordinates::Vector{<:Coordinate{RealT}}, groupid_vec::Vector{<:Integer},
    box_size::Coordinate{<:Real};
    x::Bool = true, y::Bool = true, z::Bool = true
    )::Vector{<:Coordinate{RealT}} where RealT <: Real
    half_box_size = box_size * 0.5
    fix_pbc(coordinates, groupid_vec, box_size, half_box_size, x=x, y=y, z=z)
end

"""
    fix_pbc(trj::Trajectory, box_size::Coordinate)::Trajectory

Fix residues splited by periodic boundary condition. This is more specific version of fix_pbc function for trajectory handling.
If the residue over the boundary of the box, the atoms in the residue which separated from first atom of that move to the position where the atom locate without periodic boundary condition.
"""
function fix_pbc(trj::Trajectory{RealT}, box_size::Coordinate{<:Real};
    x::Bool = true, y::Bool = true, z::Bool = true
    )::Trajectory{RealT} where RealT <: Real

    new_trj = deepcopy(trj)
    half_box_size = box_size * 0.5

    resid_vec  = map(attr->attr.resid, new_trj.attributes)
    unique_resid_vec = unique(resid_vec)
    coordinates = new_trj.coordinates
    attributes  = new_trj.attributes
    for resid in unique_resid_vec
        same_mol_indices = findall(id->id==resid, resid_vec)
        @views sbj_coords = coordinates[same_mol_indices[2:end], :]
        @views dist2first_mat = sbj_coords .- permutedims(coordinates[same_mol_indices[1], :])
        if x
            for (coord, dist) in zip(sbj_coords, dist2first_mat)
                coord.x = abs(dist.x) < half_box_size.x ? coord.x : coord.x - sign(dist.x) * box_size.x
            end
        end

        if y
            for (coord, dist) in zip(sbj_coords, dist2first_mat)
                coord.y = abs(dist.y) < half_box_size.y ? coord.y : coord.y - sign(dist.y) * box_size.y
            end
        end

        if z
            for (coord, dist) in zip(sbj_coords, dist2first_mat)
                coord.z = abs(dist.z) < half_box_size.z ? coord.z : coord.z - sign(dist.z) * box_size.z
            end
        end
    end
    new_trj
end

function fix_pbc!(trj::Trajectory{RealT}, box_size::Coordinate{<:Real};
    x::Bool = true, y::Bool = true, z::Bool = true
    ) where RealT <: Real

    # fix atom in over pbc residue based on first atom of the residue
    half_box_size = box_size * 0.5

    resid_vec  = map(attr->attr.resid, trj.attributes)
    unique_resid_vec = unique(resid_vec)
    coordinates = trj.coordinates
    attributes  = trj.attributes
    for resid in unique_resid_vec
        same_mol_indices = findall(id->id==resid, resid_vec)
        @views sbj_coords = coordinates[same_mol_indices[2:end], :]
        @views dist2first_mat = sbj_coords .- permutedims(coordinates[same_mol_indices[1], :])
        if x
            for (coord, dist) in zip(sbj_coords, dist2first_mat)
                coord.x = abs(dist.x) < half_box_size.x ? coord.x : coord.x - sign(dist.x) * box_size.x
            end
        end
        if y
            for (coord, dist) in zip(sbj_coords, dist2first_mat)
                coord.y = abs(dist.y) < half_box_size.y ? coord.y : coord.y - sign(dist.y) * box_size.y
            end
        end
        if z
            for (coord, dist) in zip(sbj_coords, dist2first_mat)
                coord.z = abs(dist.z) < half_box_size.z ? coord.z : coord.z - sign(dist.z) * box_size.z
            end
        end
    end
end

"""
    move_pbc_center(coordinates::Vector{Coordinate},
                    new_center::Coordinate, box_size::Coordinate)
                    )::Vector{Coordinate}
Fix coordinate to consistent with new box center and box size.
Box size should match to original periodic boundary box size.
"""
function move_pbc_center(coordinates::Vector{<:Coordinate{RealT}},
    center::Coordinate{<:Real}, box_size::Coordinate{<:Real}
    )::Vector{Coordinate{RealT}} where RealT <: Real

    new_coords = deepcopy(coordinates)
    half_box = box_size * 0.5
    for coord in new_coords
        dist = coord - center
        coord.x = abs(dist.x) < half_box.x ? coord.x : coord.x - sign(dist.x) * box_size.x
        coord.y = abs(dist.y) < half_box.y ? coord.y : coord.y - sign(dist.y) * box_size.y
        coord.z = abs(dist.z) < half_box.z ? coord.z : coord.z - sign(dist.z) * box_size.z
    end
    new_coords
end
