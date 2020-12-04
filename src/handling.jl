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
                            upper_bound::Coordinate, lower_bound::Coordinate)
    ::Matrix{Real}

Calculate distance matrix for all combination between `coords1` and `coord2` considering periodic boundary condition.
The box size is specified by two edges, `lower_bound` and `upper_bound`, each corresponds to the lower left front and the upper right back.
The row of returned matrix coresspond to `coords1`, and the column correspond to `coords2`.
"""
function pair_length_matrix_pbc(first_coords::ArrayT1, second_coords::ArrayT2,
    upper_bound::Coordinate{<:Real}, lower_bound::Coordinate{<:Real}
    )::Matrix{<:Real} where {ArrayT1 <: AbstractArray{<:Coordinate{<:Real}, 1},
                             ArrayT2 <: AbstractArray{<:Coordinate{<:Real}, 1}}

    distance_pbc.(first_coords, reshape(second_coords, (1, length(second_coords))),
                  upper_bound, lower_bound)
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
    pair_length_matrix_pbc(trj::Trajectory,
                           upper_bound::Coordinate, lower_bound::Coordinate;
                           frame_indices::Union{Vector, OrdinalRange, Colon}       = :,
                           first_atom_indices::Union{Vector, OrdinalRange, Colon}  = :,
                           second_atom_indices::Union{Vector, OrdinalRange, Colon} = atom_indices1)
    ::Vector{Matrix{Coordinate}}

Calculate distance matrix for all combinations between `first_atom_indices` and `second_atom_indices` considering periodic boundary condition.
The box size is specified by two edges, `lower_bound` and `upper_bound`, each corresponds to the lower left front and the upper right back.
This cauculation apply to each frame of trajectory and the result matrices are stored in Vector.
The target frame can be restricted by pass indeces vector or range to `frame_indices`.
"""
function pair_length_matrix_pbc(trj::Trajectory{<:Real},
    upper_bound::Coordinate{<:Real}, lower_bound::Coordinate{<:Real};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Vector{<:Matrix{<:Real}}

    zip_iterate4frame = zip(eachcol(view(trj.coordinates, first_atom_indices, frame_indices)),
                            eachcol(view(trj.coordinates, second_atom_indices, frame_indices)))
    map(coords_vec_pair -> pair_length_matrix_pbc(coords_vec_pair...,
                                                  upper_bound, lower_bound),
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
    contact_bool_matrix_pbc(threshold::Real, trj::Trajectory
                            upper_bound::Coordinate, lower_bound::Coordinate;
                            frame_indices::Union{Array, OrdinalRange, Colon}       = :,
                            first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,
                            second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)
    ::Vector{Matrix{Bool}}

Judge contact is formed or not considering periodic boundary condition. If the distance between two coordinate is shorter than threshold, contact is considered to be formed. In returned vector of matrices, each matrix correspond to contact matrix of each frame.
The box size is specified by two edges, `lower_bound` and `upper_bound`, each corresponds to the lower left front and the upper right back.
You can specify the target frames or atoms by `frame_indices`, `first_atom_indices` or `second_atom_indices`. When you specify the target atoms, the row of matrices corresponds to first_atom_indices and column of matrices corresponds to second_atom_indices.
"""

function contact_bool_matrix_pbc(threshold::RealT, trj::Trajectory{<:Real},
    upper_bound::Coordinate{<:Real}, lower_bound::Coordinate{<:Real};
    frame_indices::Union{Vector, OrdinalRange, Colon} = :,
    first_atom_indices::Union{Vector, OrdinalRange, Colon} = :,
    second_atom_indices::Union{Vector, OrdinalRange, Colon} = first_atom_indices
    )::Vector{Matrix{Bool}} where RealT <: Real

    length_mat_arr = pair_length_matrix_pbc(trj, upper_bound, lower_bound,
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
                 upper_bound::Coordinate, lower_bound::Coordinate)
    ::Real

Calculate distance `first_atom` and `second_atom` considering periodic boundary condition.
The box size is specified by two edges, `lower_bound` and `upper_bound`, each corresponds to the lower left front and the upper right back.
"""
function distance_pbc(first_coord::Coordinate{RealT},  second_coord::Coordinate{RealT},
    upper_bound::Coordinate{<:Real}, lower_bound::Coordinate{<:Real}
    )::RealT where RealT <: Real

    box_vec  = upper_bound - lower_bound
    dist_vec = first_coord  - second_coord
    x = abs(dist_vec.x) < box_vec.x * 0.5 ? dist_vec.x : box_vec.x - abs(dist_vec.x)
    y = abs(dist_vec.y) < box_vec.y * 0.5 ? dist_vec.y : box_vec.y - abs(dist_vec.y)
    z = abs(dist_vec.z) < box_vec.z * 0.5 ? dist_vec.z : box_vec.z - abs(dist_vec.z)
    sqrt(x^2 + y^2 + z^2)
end

function fix_pbc!(trj::Trajectory{RealT},
    lower_bound::Coordinate{RealT}, upper_bound::Coordinate{RealT}
    ) where RealT <: Real

    # TODO: remove side effect
    # fix atom in over pbc residue based on first atom of the residue
    box_coord = upper_bound - lower_bound
    half_box_coord = box_coord * 0.5

    resid_vec  = map(attr->attr.resid, trj.attributes)
    unique_resid_vec = unique(resid_vec)
    coordinates = trj.coordinates
    attributes  = trj.attributes
    for resid in unique_resid_vec
        same_mol_indices = findall(id->id==resid, resid_vec)
        @views sbj_coords = coordinates[same_mol_indices[2:end], :]
        @views dist2first_mat = sbj_coords .- permutedims(coordinates[same_mol_indices[1], :])
        for (coord, dist) in zip(sbj_coords, dist2first_mat)
            coord.x = abs(dist.x) < half_box_coord.x ? coord.x : coord.x - sign(dist.x) * box_coord.x
            coord.y = abs(dist.y) < half_box_coord.y ? coord.y : coord.y - sign(dist.y) * box_coord.y
            coord.z = abs(dist.z) < half_box_coord.z ? coord.z : coord.z - sign(dist.z) * box_coord.z
        end
    end
end
