module MolHandler

export Coordinate, norm, distance
export Atom, Attribute, Frame, Trajectory
export read_dcd, write_dcd, read_pdb, write_pdb, read_xyz, write_xyz
export get_frame, get_atom, clip_trajectory, center_of_mass,
       pair_length_matrix, pair_length_matrix_pbc, pair_length_matrix_parallel,
       contact_bool_matrix, contact_bool_matrix_pbc, contact_bool_matrix_parallel,
       contact_probability_matrix, contact_probability_matrix_parallel,
       radius_of_gyration,
       atom_mass, residue_mass,
       distance_pbc, fix_pbc, fix_pbc!, move_pbc_center

# codes
include("variables.jl")
include("coordinate.jl")
include("component.jl")
include("fileio.jl")
include("handling.jl")

end # module
