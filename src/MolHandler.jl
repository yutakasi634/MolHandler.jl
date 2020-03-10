module MolHandler

export Coordinate, norm, distance
export Atom, Attribute, Frame, Trajectory
export readdcd, readpdb
export get_frame, get_atom, clip_trajectory, center_of_mass, pair_length_matrix, contact_bool_matrix, contact_probability_matrix

# codes
include("coordinate.jl")
include("component.jl")
include("fileio.jl")
include("handling.jl")

end # module
