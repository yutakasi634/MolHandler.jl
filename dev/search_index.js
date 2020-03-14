var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions-1","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/#Index-1","page":"Functions","title":"Index","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"Order = [:function]","category":"page"},{"location":"functions/#File-IO-1","page":"Functions","title":"File IO","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"read_dcd\nwrite_dcd\nread_pdb","category":"page"},{"location":"functions/#MolHandler.read_dcd","page":"Functions","title":"MolHandler.read_dcd","text":"read_dcd(filename::String)::Trajectory\n\nReturn Trajectory object which filled coordinates, nframe, natom fields.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.write_dcd","page":"Functions","title":"MolHandler.write_dcd","text":"write_dcd(filename::String, trj::Trajectory;\n          save_step::Integer = 1, total_step::Integer = trj.nframe,\n          unit_num::Integer  = 1,  time_step::Real = 1.0f0)\n\nWrite coordinates, nframe and natom infomation of trj to dcd file, named filename. dcd file contain other information - save interval step, unit(chain) number in the system, total step and time step of original trajectory -, and if you do not specify these parameter, these are filled with default value.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.read_pdb","page":"Functions","title":"MolHandler.read_pdb","text":"read_pdb(filename::String)::Trajectory\n\nReturn Trajectory object which filled all field.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Trajectory-handling-1","page":"Functions","title":"Trajectory handling","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"get_frame\nget_atom\nclip_trajectory\ncenter_of_mass\npair_length_matrix\ncontact_bool_matrix\ncontact_probability_matrix","category":"page"},{"location":"functions/#MolHandler.get_frame","page":"Functions","title":"MolHandler.get_frame","text":"get_frame(frame_idx::Int64, trajectory::Trajectory)\n::Frame\n\nReturn Frame object which correspond to framd_idx frame from trajectory.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.get_atom","page":"Functions","title":"MolHandler.get_atom","text":"get_atom(query::Int64, trajectory::Trajectory)\n::Vector{Atom}\n\nReturn Vector of queryth Atom object from trajectory.\n\n\n\n\n\nget_atom(query::Int64, frame::Frame)\n::Atom\n\nReturn queryth Atom object from frame.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.clip_trajectory","page":"Functions","title":"MolHandler.clip_trajectory","text":"clip_trajectory(query::Union{Integer, Array{T, 1}, OrdinalRange}, trj::Trajectory;\n                query_key::Symbol = :frame)\n::Trajecotory where T <: Integer\n\nClip a part of trajectory you specified. Expected query key is one of the following\n\n:frame : in this case, query specify the frame index.\n:atom  : in this case, query specify the atom index.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.center_of_mass","page":"Functions","title":"MolHandler.center_of_mass","text":"center_of_mass(query::Trajectory;\n               indices::Union{Array, OrdinalRange, Colon} = :,\n               geometric::Bool = false)\n::Vector{Coordinate{Real}}\n\nCalculate the center of mass of trajectory for specified atom indices. If you set geometric is true, this function calculate geometric center of mass.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.pair_length_matrix","page":"Functions","title":"MolHandler.pair_length_matrix","text":"pair_length_matrix(corrds1::Vector{Coordinate}, coords2::Vector{Coordinate})\n::Matrix{Coordinate}\n\nCalculate distance matrix for all combination between coords1 and coord2. The row of returned matrix coresspond to coords1, and the column correspond to coords2.\n\n\n\n\n\npair_length_matrix(trj::Trajectory;\n                   frame_indices::Union{Array, OrdinalRange, Colon}       = :,\n                   first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,\n                   second_atom_indices::Union{Array, OrdinalRange, Colon} = atom_indices1)\n::Vector{Matrix{Coordinate}}\n\nCalculate distance matrix for all combination between first_atom_indices and second_atom_indices. This cauculation apply to each frame of trajectory and the result matrices are stored in Vector. The target frame can be restricted by pass indeces vector or range to frame_indices.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.contact_bool_matrix","page":"Functions","title":"MolHandler.contact_bool_matrix","text":"contact_bool_matrix(threshold::Real, trj::Trajectory;\n                    frame_indices::Union{Array, OrdinalRange, Colon}       = :,\n                    first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,\n                    second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)\n::Vector{Matrix{Bool}}\n\nJudge contact is formed or not. If the distance between two coordinate is shorter than threshold, contact is considered to be formed. In returned matrix, the row of matrices corresponds to firstatomindices and column of matrices corresponds to secondatomindices.\n\n\n\n\n\ncontact_bool_matrix(threshold::Real, trj::Trajectory;\n                    frame_indices::Union{Array, OrdinalRange, Colon}       = :,\n                    first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,\n                    second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)\n::Vector{Matrix{Bool}}\n\nJudge contact is formed or not. If the distance between two coordinate is shorter than threshold, contact is considered to be formed. In returned vector of matrices, each matrix correspond to contact matrix of each frame. You can specify the target frames or atoms by frame_indices, first_atom_indices or second_atom_indices. When you specify the target atoms, the row of matrices corresponds to firstatomindices and column of matrices corresponds to secondatomindices.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.contact_probability_matrix","page":"Functions","title":"MolHandler.contact_probability_matrix","text":"contact_probability_matrix(threshold::Real, trj::Trajectory;\n                           frame_indices::Union{Array, OrdinalRange, Colon}       = :,\n                           first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,\n                           second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)\n::Matrix{Real}\n\nCalculate contact formation probability over the trajectory. If the distance between two coordinate is shorter than threshold, contact is considered to be formed. You can specify the target frames or atoms by frame_indices, first_atom_indices or second_atom_indices. When you specify the target atoms, the row of matrices corresponds to firstatomindices and column of matrices corresponds to secondatomindices.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Multi-Threading-1","page":"Functions","title":"Multi-Threading","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"pair_length_matrix_parallel\ncontact_bool_matrix_parallel\ncontact_probability_matrix_parallel","category":"page"},{"location":"functions/#MolHandler.pair_length_matrix_parallel","page":"Functions","title":"MolHandler.pair_length_matrix_parallel","text":"pair_length_matrix_parallel(trj::Trajectory;\n                            frame_indices::Union{Array, OrdinalRange, Colon}       = :,\n                            first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,\n                            second_atom_indices::Union{Array, OrdinalRange, Colon} = atom_indices1)\n::Vector{Matrix{Coordinate}}\n\nMulti-threads version of pairlengthmatrix. If you set available threads number to Threads.nthreads(), this function would faster than non-parallel version.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.contact_bool_matrix_parallel","page":"Functions","title":"MolHandler.contact_bool_matrix_parallel","text":"contact_bool_matrix_parallel(threshold::Real, trj::Trajectory;\n                             frame_indices::Union{Array, OrdinalRange, Colon}       = :,\n                             first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,\n                             second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)\n::Vector{Matrix{Bool}}\n\nMulti-threads version of contactboolmatrix. If you set available threads number to Threads.nthreads(), this function would faster than non-parallel version.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.contact_probability_matrix_parallel","page":"Functions","title":"MolHandler.contact_probability_matrix_parallel","text":"contact_probability_matrix_parallel(threshold::Real, trj::Trajectory;\n                                    frame_indices::Union{Array, OrdinalRange, Colon}       = :,\n                                    first_atom_indices::Union{Array, OrdinalRange, Colon}  = :,\n                                    second_atom_indices::Union{Array, OrdinalRange, Colon} = first_atom_indices)\n::Matrix{Real}\n\nMulti-threads version of contactprobabilitymatrix. If you set available threads number to Threads.nthreads(), this function would faster than non-parallel version.\n\n\n\n\n\n","category":"function"},{"location":"structs/#Structs-1","page":"Structs","title":"Structs","text":"","category":"section"},{"location":"structs/#Index-1","page":"Structs","title":"Index","text":"","category":"section"},{"location":"structs/#","page":"Structs","title":"Structs","text":"Order = [:type]","category":"page"},{"location":"structs/#","page":"Structs","title":"Structs","text":"Trajectory\nCoordinate\nAttribute\nFrame\nAtom","category":"page"},{"location":"structs/#MolHandler.Trajectory","page":"Structs","title":"MolHandler.Trajectory","text":"This struct corresponding to MD trajectory. This contains 4 fields like below.\n\ncoordinates : Matrix{Coordinate}\nattributes  : Vector{Attribute}\nnatom       : Number of atoms.\nnframe      : Number of frames in the trajectory.\n\nThe column of coordinates matrix means one snapshot. The row of coordinates matrix means time series of one atom.\n\n[ a d g j\n  b e h k\n  c f i l ]\n\nIn this case, one snapshot correspond to [a b c].\n\n\n\n\n\n","category":"type"},{"location":"structs/#MolHandler.Coordinate","page":"Structs","title":"MolHandler.Coordinate","text":"This correspond to coordinate of atom. This contain 3 fields like below.\n\nx <: Real\ny <: Real\nz <: Real\n\nSome operators and functions were overloaded to make Coordinate objects broadcastable.\n\n\n\n\n\n","category":"type"},{"location":"structs/#MolHandler.Attribute","page":"Structs","title":"MolHandler.Attribute","text":"This struct mean all atomic information other than coordinates. This contains 5 fields like below.\n\nresname  : Union{String,  Nothing}\nresid    : Union{Int64,   Nothing}\natomname : Union{String,  Nothing}\natomid   : Union{Int64,   Nothing}\nmass     : Union{Float32, Nothing}\n\n\n\n\n\n","category":"type"},{"location":"structs/#MolHandler.Frame","page":"Structs","title":"MolHandler.Frame","text":"This struct correspond to one snapshot of trajectory. This contains 3 fields like below.\n\ncoordinates : Vector{Coordinate}\nattribute   : Vector{Attribute}\nnatom       : Number of atoms.\n\n\n\n\n\n","category":"type"},{"location":"structs/#MolHandler.Atom","page":"Structs","title":"MolHandler.Atom","text":"This struct have all information of corresponding atom. This contains 2 fields like below.\n\ncoordinate : Reference to Coordinate object.\nattribute  : Reference to Attribute object.\n\n\n\n\n\n","category":"type"},{"location":"#MolHandler.jl-Documentation-1","page":"Top","title":"MolHandler.jl Documentation","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"This package is for handling MD trajectory data. Supported format is dcd and pdb.","category":"page"},{"location":"#How-to-install-1","page":"Top","title":"How to install","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"This package assumes julia 1.2.","category":"page"},{"location":"#","page":"Top","title":"Top","text":"pkg> add https://github.com/yutakasi634/MolHandler.jl.git\njulia> using MolHandler","category":"page"},{"location":"#Example-of-use-1","page":"Top","title":"Example of use","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"julia> using MolHandler\njulia> trj = read_dcd(\"trajectory.dcd\")\njulia> #trj = read_pdb(\"structure.pdb\")\njulia> trj.coordinates[:,1] # get first snapshot as Coordinate object array.\njulia> trj.coordinates[1,:] # get first atom coordinate time series by Coordinate object array.\njulia> frame      = get_frame(1, trj) # get first frame as Frame object.\njulia> atom_array = get_atom(1, trj) # get first atom time series as Atom array.\njulia> atom       = get_atom(2, frame) # get second atom as Atom object.","category":"page"},{"location":"#[Functions](@ref)-list-1","page":"Top","title":"Functions list","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"Order = [:function]","category":"page"},{"location":"#[Structs](@ref)-list-1","page":"Top","title":"Structs list","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"Order = [:type]","category":"page"}]
}
