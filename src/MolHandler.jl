module MolHandler

export Coordinate, norm, distance
export Atom, Attribute, Frame, Trajectory
export read_dcd, write_dcd, read_pdb
export get_frame, get_atom, clip_trajectory, center_of_mass, pair_length_matrix, pair_length_matrix_parallel, contact_bool_matrix, contact_bool_matrix_parallel, contact_probability_matrix, contact_probability_matrix_parallel

# codes
include("coordinate.jl")
include("component.jl")
include("fileio.jl")
include("handling.jl")

end # module
