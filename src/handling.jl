import Base.Threads
import LinearAlgebra
import Random

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
        Trajectory(reshape(trj.coordinates[:, query], (trj.natom, 1)), attributes=trj.attributes)
    elseif query_key == :atom
        Trajectory(reshape(trj.coordinates[query, :], (1, trj.nframe)),
                   attributes=[trj.attributes[query]])
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
        Trajectory(trj.coordinates[:, query], attributes=trj.attributes)
    elseif query_key == :atom
        Trajectory(trj.coordinates[query, :], attributes=trj.attributes[query])
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
    center_of_mass(coordinates::Vector{Coordiante}, mass_vec{Real};
                   atom_indices::Union{Vector, OrdinalRange, Colon} = :)
    ::Coordinate{RealT}

Calculate the (non-geometric) center of mass of Coordinate vector for specified atom indices.
You have to pass mass information vector which have same length to coordinates vector.
"""
function center_of_mass(coordinates::Vector{<:Coordinate{RealT}}, mass_vec::Vector{<:Real};
    atom_indices::Union{Vector, OrdinalRange, Colon} = :)::Coordinate{RealT} where RealT <: Real

    if length(coordinates) != length(mass_vec)
        throw(AssertionError("""
                             coordinate vector and mass vector should be same length.
                             """))
    end

    selected_coord = view(coordinates, atom_indices)
    selected_mass  = view(mass_vec,    atom_indices)
    sum(selected_coord .* selected_mass) / sum(selected_mass)
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
    box_size::Coordinate{<:Real};
    x::Bool = true, y::Bool = true, z::Bool = true)::Vector{Coordinate}

Fix the atom group splited by periodic boundary box. This group is defined by `groupid_vec` which contain each particles group id.
For example, if your system have 3 atom, and atom 1 and 2 are group 1, and atom 3 is group 2, this `groupid_vec` is [1, 1, 2].
If the group over the boundary of the box, the atoms in the group which separated from first atom of that move to the position where the atom locate without periodic boundary condition.
You can specify the axies which applied fixing by x,y and option. If you set `z = false`, only x and y coordinate fixed and z don't change.
"""
function fix_pbc(coordinates::Vector{<:Coordinate{RealT}}, groupid_vec::Vector{<:Integer},
    box_size::Coordinate{<:Real};
    x::Bool = true, y::Bool = true, z::Bool = true
    )::Vector{<:Coordinate{RealT}} where RealT <: Real
    half_box_size = box_size * 0.5
    fix_pbc(coordinates, groupid_vec, box_size, half_box_size, x=x, y=y, z=z)
end

# This function get half_box_size, because if process will apply to all trajectory,
# half box size calculation will occur every time frame, this is overhead.
function fix_pbc!(coordinates::ArrayT, contact_matrix::Matrix{Bool},
    box_size::Coordinate{<:Real}, half_box_size::Coordinate{<:Real};
    x::Bool= true, y::Bool = true, z::Bool = true
    ) where ArrayT <: AbstractArray{<:Coordinate{<:Real}, 1}

    if length(coordinates) != size(contact_matrix)[1] || length(coordinates) != size(contact_matrix)[2]
        throw(AssertionError("""
                             the size of coordinate vector and that of each dimension of contact matrix should be same length.
                             """))
    end

    fixed_atom_arr = fill(false, length(coordinates))
    stack_arr      = Array{Int64, 1}()

    for atom_idx in 1:length(coordinates)
        # if atom already fixed, go to next atom
        if fixed_atom_arr[atom_idx]
            continue
        end

        fixed_atom_arr[atom_idx] = true
        push!(stack_arr, atom_idx)
        while !isempty(stack_arr)
            ref_atom_idx = popfirst!(stack_arr)
            ref_atom     = coordinates[ref_atom_idx]
            for (target_idx, is_target_atom) in enumerate(view(contact_matrix, ref_atom_idx, :))
                if !is_target_atom || fixed_atom_arr[target_idx]
                    continue
                end
                fixed_atom_arr[target_idx] = true
                push!(stack_arr, target_idx)
                target_atom = coordinates[target_idx]
                dist        = target_atom - ref_atom

                if x
                    target_atom.x = abs(dist.x) < half_box_size.x ? target_atom.x : target_atom.x - sign(dist.x) * box_size.x
                end

                if y
                    target_atom.y = abs(dist.y) < half_box_size.y ? target_atom.y : target_atom.y - sign(dist.y) * box_size.y
                end

                if z
                    target_atom.z = abs(dist.z) < half_box_size.z ? target_atom.z : target_atom.z - sign(dist.z) * box_size.z
                end
            end
        end
    end
end


"""
    fix_pbc!(coordinates::Vector{Coordinate}, contact_matrix::Matrix{Bool},
    box_size::Coordinate{<:Real};
    x::Bool = true, y::Bool = true, z::Bool = true)::Vector{Coordinate}

Fix the atom group splited by periodic boundary box. This group is juged by contact_matrix.
If the group over the boundary of the box, the atoms in the group which separated from first atom of that move to the position where the atom locate without periodic boundary condition.
You can specify the axies which applied fixing by x,y and option. If you set `z = false`, only x and y coordinate fixed and z don't change.
"""
function fix_pbc!(coordinates::ArrayT, contact_matrix::Matrix{Bool},
    box_size::Coordinate{<:Real};
    x::Bool = true, y::Bool = true, z::Bool = true
    ) where ArrayT <: AbstractArray{<:Coordinate{<:Real}, 1}
    half_box_size = box_size * 0.5
    fix_pbc!(coordinates, contact_matrix, box_size, half_box_size, x=x, y=y, z=z)
end

"""
    fix_pbc(trj::Trajectory, groupid_vec::Vector{Integer},
    box_size::Coordinate;
    x::Bool = true, y::Bool = true, z::Bool = true)::Trajectory

Fix residues splited by periodic boundary condition. This is more specific version of `fix_pbc` function for trajectory handling.
For example, if your system have 3 atom, and atom 1 and 2 are group 1, and atom 3 is group 2, this `groupid_vec` is [1, 1, 2].
If the group over the boundary of the box, the atoms in the group which more than half a box away from the first atom of that move to the position where the atom locate without periodic boundary condition.
You can specify the axies which applied fixing by x,y and option. If you set `z = false`, only x and y coordinate fixed and z don't change.
"""
function fix_pbc(trj::Trajectory{RealT}, groupid_vec::Vector{<:Integer},
    box_size::Coordinate{<:Real};
    x::Bool = true, y::Bool = true, z::Bool = true
    )::Trajectory{RealT} where RealT <: Real

    if size(trj.coordinates)[1] != length(groupid_vec)
        throw(AssertionError("""
                             groupid_vec do not have appropriate length. This length is $(length(groupid_vec)), but the number of particle in trajectory is $(size(trj.coordinates)[1]).
                             """))
    end

    new_trj = deepcopy(trj)
    half_box_size = box_size * 0.5

    unique_groupid_vec = unique(groupid_vec)
    coordinates        = new_trj.coordinates
    for groupid in unique_groupid_vec
        same_group_indices = findall(id->id==groupid, groupid_vec)
        @views sbj_coords = coordinates[same_group_indices[2:end], :]
        @views dist2first_mat = sbj_coords .- permutedims(coordinates[same_group_indices[1], :])
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

function fix_pbc(trj::Trajectory{RealT}, groupid_vec::Vector{<:Integer};
    x::Bool = true, y::Bool = true, z::Bool = true
    )::Trajectory{RealT} where RealT <: Real

    if size(trj.coordinates)[2] != length(trj.boxes)
        throw(AssertionError("""
                             Trajctory object do not have appropriate box information. The number of frame is $(size(trj.coordinates)[2]), but the number of box information is $(length(trj.boxes)).
                             """))
    end

    unique_groupid_vec = unique(groupid_vec)

    new_trj = deepcopy(trj)
    coordinates = new_trj.coordinates
    boxes       = new_trj.boxes

    for frame_idx in 1:size(coordinates)[2]
        box_size      = boxes[frame_idx]
        half_box_size = box_size * 0.5
        for groupid in unique_groupid_vec
            same_group_indices = findall(id->id==groupid, groupid_vec)
            @views sbj_coords = coordinates[same_group_indices[2:end], frame_idx]
            @views dist2first = sbj_coords .- coordinates[same_group_indices[1], frame_idx]
            if x
                for (coord, dist) in zip(sbj_coords, dist2first)
                    coord.x = abs(dist.x) < half_box_size.x ? coord.x : coord.x - sign(dist.x) * box_size.x
                end
            end
            if y
                for (coord, dist) in zip(sbj_coords, dist2first)
                    coord.y = abs(dist.y) < half_box_size.y ? coord.y : coord.y - sign(dist.y) * box_size.y
                end
            end
            if z
                for (coord, dist) in zip(sbj_coords, dist2first)
                    coord.z = abs(dist.z) < half_box_size.z ? coord.z : coord.z - sign(dist.z) * box_size.z
                end
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

function fix_pbc_along_time(trj::Trajectory{RealT};
    x::Bool = true, y::Bool = true, z::Bool = true
    )::Trajectory{RealT} where RealT <: Real

    new_trj     = deepcopy(trj)
    coordinates = new_trj.coordinates
    boxes       = new_trj.boxes

    @views previous_coords = coordinates[:, 1]
    for frame_idx in 2:size(coordinates)[2]
        box_size      = boxes[frame_idx]
        half_box_size = box_size * 0.5

        @views sbj_coords = coordinates[:, frame_idx]
        dist2previous = sbj_coords .- previous_coords

        if x
            for (coord, dist) in zip(sbj_coords, dist2previous)
                if abs(dist.x) > half_box_size.x
                    fix_coef =
                        sign(dist.x) *
                        (div((abs(dist.x) - half_box_size.x), box_size.x) + 1)
                    coord.x = coord.x - fix_coef * box_size.x
                end
            end
        end
        if y
            for (coord, dist) in zip(sbj_coords, dist2previous)
                if abs(dist.y) > half_box_size.y
                    fiy_coef =
                        sign(dist.y) *
                        (div((abs(dist.y) - half_box_size.y), box_size.y) + 1)
                    coord.y = coord.y - fiy_coef * box_size.y
                end
            end
        end
        if z
            for (coord, dist) in zip(sbj_coords, dist2previous)
                if abs(dist.z) > half_box_size.z
                    fiz_coef =
                        sign(dist.z) *
                        (div((abs(dist.z) - half_box_size.z), box_size.z) + 1)
                    coord.z = coord.z - fiz_coef * box_size.z
                end
            end
        end

        @views previous_coords = sbj_coords
    end
    new_trj
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

"""
    sasa(coordinates::Vector{Coordinate}, radiuses::Vector{Real},
         solvent_radius::RealT = 1.4,     MCtrial_num::Integer = 1000
         seed::Integer         = undef)::Vector{Real}
Calculate Solvent Accessible Surface Area(SASA) for each particle.
The length of radiuses should match to the number of particles in coordinates.
The default value of solvent radius, 1.4, corresponds to water case 1.4 Å.
MCtrial means MonteCarlo trial num for calculate the solvent accessible area of each particle.
For each particle, this function generate the point MCtrial times on the solvent accessible surface and count the number of non-overlapping dots with other particle.
"""
function sasa(coordinates::Vector{Coordinate{RealT}}, radiuses::Vector{RealT},
    solvent_radius::RealT = RealT(1.4), mctrial_num::Integer = 1000,
    seed                  = undef
    )::Vector{RealT} where RealT <: Real

    @assert length(coordinates) == length(radiuses) "the length of coordinates vector and radiuses vector must be same."

    rng = Random.MersenneTwister()
    if seed != undef
        rng.seed!(seed)
    end

    sasa_arr   = Vector{RealT}()
    n_particle = length(coordinates)
    for sbj_coord_idx in 1:n_particle
        # search neigbour particles
        sbj_coord = coordinates[sbj_coord_idx]
        sbj_radi  = radiuses[sbj_coord_idx]
        neibour_indices = Vector{Integer}()
        for coord_idx in 1:n_particle
            coord = coordinates[coord_idx]
            dist  = norm(sbj_coord - coord)
            if dist < sbj_radi + radiuses[coord_idx] + solvent_radius*2
                push!(neibour_indices, coord_idx)
            end
        end

        access_count = 0
        sa_radi      = sbj_radi + solvent_radius
        for traial_idx in 1:mctrial_num
            u   = Random.rand(rng)
            phi = 2.0*π*Random.rand(rng)
            sin_theta = 1.0 - u^2

            dot = Coordinate(sa_radi*sin_theta*cos(phi),
                             sa_radi*sin_theta*sin(phi),
                             sa_radi*u) + sbj_coord
            overlap = false
            for coord_idx in neibour_indices
                dist = norm(dot - coordinates[coord_idx])
                if dist < radiuses[coord_idx] + solvent_radius
                    overlap = true
                    break
                end
            end
            if overlap
                access_count += 1
            end
        end
        sa_ratio = access_count / mctrial_num
        sasa     = 4.0*π*sa_radi^2 * sa_ratio
        push!(sasa_arr, sasa)
    end
    sasa_arr
end

"""
    rotate(coordinates::Coordinate,
           rotate_x::Real, rotate_y::Real, rotate_z::Real
           )::Vector{Real}
Rotate particles coordinate around x with rotate_x, y with rotate_y, z with rotate_z. The unit is radian.
"""
function rotate(coordinate::Coordinate{RealT},
    rotate_x::RealT, rotate_y::RealT, rotate_z::RealT,
    )::Coordinate{RealT} where RealT <: Real
    x = coordinate.x
    y = coordinate.y
    z = coordinate.z
    sin_x = sin(rotate_x)
    sin_y = sin(rotate_y)
    sin_z = sin(rotate_z)
    cos_x = cos(rotate_x)
    cos_y = cos(rotate_y)
    cos_z = cos(rotate_z)
    new_x = x*cos_z*cos_y + y*(cos_z*sin_y*cos_x - sin_z*cos_x) + z*(cos_z*sin_y*cos_x + sin_z*sin_x)
    new_y = x*sin_z*cos_y + y*(sin_z*sin_y*sin_x + cos_z*cos_x) + z*(sin_z*sin_y*cos_x - cos_z*sin_x)
    new_z = -x*sin_y + y*cos_y*sin_x + z*cos_y*cos_x

    Coordinate(new_x, new_y, new_z)
end

"""
    rotate(trj::Trajectory{Real},
           rotate_x::Real, rotate_y::Real, rotate_z::Real
           )::Trajectory{Real}
Rotate all particles coordinate in trajectory around x with rotate_x, y with rotate_y, z with rotate_z. The unit is radian.
"""
function rotate(trj::Trajectory{RealT},
    rotate_x::RealT, rotate_y::RealT, rotate_z::RealT
    )::Trajectory{RealT} where RealT <: Real

    new_trj     = deepcopy(trj)
    coordinates = new_trj.coordinates

    for idx in eachindex(coordinates)
        coordinates[idx] = rotate(coordinates[idx],
                                  rotate_x, rotate_y, rotate_z)
    end
    new_trj
end
